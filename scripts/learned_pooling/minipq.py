"""Minimal single-column Parquet reader (BYTE_ARRAY, dictionary or plain).
Dependency-free: thrift compact protocol + snappy/gzip + RLE hybrid decoding."""
import struct, zlib

def _uvarint(b, i):
    shift = 0; out = 0
    while True:
        c = b[i]; i += 1
        out |= (c & 0x7F) << shift
        if not (c & 0x80):
            return out, i
        shift += 7

def _svarint(b, i):
    v, i = _uvarint(b, i)
    return (v >> 1) ^ -(v & 1), i

def read_struct(b, i):
    out = {}; last = 0
    while True:
        h = b[i]; i += 1
        if h == 0:
            return out, i
        delta = (h & 0xF0) >> 4; ftype = h & 0x0F
        if delta == 0:
            fid, i = _svarint(b, i)
        else:
            fid = last + delta
        last = fid
        val, i = read_value(b, i, ftype)
        out[fid] = val

def read_value(b, i, ftype):
    if ftype == 1: return True, i
    if ftype == 2: return False, i
    if ftype == 3:
        v = b[i]; return (v - 256 if v > 127 else v), i + 1
    if ftype in (4, 5, 6): return _svarint(b, i)
    if ftype == 7: return struct.unpack_from("<d", b, i)[0], i + 8
    if ftype == 8:
        n, i = _uvarint(b, i); return b[i:i + n], i + n
    if ftype in (9, 10):
        h = b[i]; i += 1
        size = (h & 0xF0) >> 4; et = h & 0x0F
        if size == 15:
            size, i = _uvarint(b, i)
        vals = []
        for _ in range(size):
            v, i = read_value(b, i, et); vals.append(v)
        return vals, i
    if ftype == 12: return read_struct(b, i)
    raise NotImplementedError("thrift type %d" % ftype)

def snappy_decompress(data):
    n, pos = _uvarint(data, 0)
    out = bytearray()
    while pos < len(data):
        tag = data[pos]; pos += 1
        t = tag & 0x03
        if t == 0:
            ln = tag >> 2
            if ln < 60:
                length = ln + 1
            else:
                nb = ln - 59
                length = int.from_bytes(data[pos:pos + nb], "little") + 1; pos += nb
            out += data[pos:pos + length]; pos += length
        else:
            if t == 1:
                length = 4 + ((tag >> 2) & 0x07)
                offset = ((tag >> 5) << 8) | data[pos]; pos += 1
            elif t == 2:
                length = (tag >> 2) + 1
                offset = int.from_bytes(data[pos:pos + 2], "little"); pos += 2
            else:
                length = (tag >> 2) + 1
                offset = int.from_bytes(data[pos:pos + 4], "little"); pos += 4
            start = len(out) - offset
            for k in range(length):
                out.append(out[start + k])
    assert len(out) == n
    return bytes(out)

def decompress(codec, data):
    if codec == 0: return data
    if codec == 1: return snappy_decompress(data)
    if codec == 2: return zlib.decompress(data, 16 + zlib.MAX_WBITS)
    raise NotImplementedError("codec %d" % codec)

def rle_hybrid(buf, bit_width, count):
    vals = []; i = 0
    if bit_width == 0:
        return [0] * count
    while len(vals) < count and i < len(buf):
        header, i = _uvarint(buf, i)
        if header & 1:
            groups = header >> 1
            nb = groups * bit_width
            bits = int.from_bytes(buf[i:i + nb], "little"); i += nb
            mask = (1 << bit_width) - 1
            for k in range(groups * 8):
                vals.append((bits >> (k * bit_width)) & mask)
        else:
            run = header >> 1
            nb = (bit_width + 7) // 8
            v = int.from_bytes(buf[i:i + nb], "little"); i += nb
            vals.extend([v] * run)
    return vals[:count]

def plain_byte_arrays(buf, count=None):
    out = []; i = 0
    while i < len(buf) and (count is None or len(out) < count):
        (n,) = struct.unpack_from("<I", buf, i); i += 4
        out.append(buf[i:i + n]); i += n
    return out

def read_column(path, column_name):
    raw = open(path, "rb").read()
    assert raw[:4] == b"PAR1" and raw[-4:] == b"PAR1"
    (flen,) = struct.unpack_from("<I", raw, len(raw) - 8)
    meta, _ = read_struct(raw, len(raw) - 8 - flen)
    values = []
    for rg in meta[4]:
        for col in rg[1]:
            md = col[3]
            if [p.decode() for p in md[3]][-1] != column_name:
                continue
            codec = md[4]; num_values = md[5]
            start = md.get(11) or md[9]
            buf = raw[start:start + md[7]]
            pos = 0; dictionary = None; got = 0
            while pos < len(buf) and got < num_values:
                ph, npos = read_struct(buf, pos)
                ptype = ph[1]; unc = ph[2]; comp = ph[3]
                page = buf[npos:npos + comp]; pos = npos + comp
                if ptype == 2:
                    dictionary = plain_byte_arrays(decompress(codec, page))
                elif ptype in (0, 3):
                    if ptype == 0:
                        d = decompress(codec, page)
                        n = ph[5][1]; enc = ph[5][2]
                        (dl_len,) = struct.unpack_from("<I", d, 0)
                        dl = rle_hybrid(d[4:4 + dl_len], 1, n)
                        body = d[4 + dl_len:]
                    else:
                        hdr = ph[8]
                        n = hdr[1]; enc = hdr[4]; dl_len = hdr[5]; rl_len = hdr[6]
                        dl = rle_hybrid(page[rl_len:rl_len + dl_len], 1, n)
                        body = decompress(codec, page[rl_len + dl_len:])
                    ndef = sum(1 for v in dl if v == 1) if dl else n
                    if enc in (2, 8):
                        vals = [dictionary[k] for k in rle_hybrid(body[1:], body[0], ndef)]
                    elif enc == 0:
                        vals = plain_byte_arrays(body, ndef)
                    else:
                        raise NotImplementedError("encoding %d" % enc)
                    it = iter(vals)
                    for v in dl:
                        values.append(next(it).decode() if v == 1 else None)
                    got += n
    return values
