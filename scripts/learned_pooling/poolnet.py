"""Learned permutation-invariant donor pooling in pure NumPy.

Three operators, all trained on donors with hand-derived gradients and Adam:
  deepsets  : mean pooling of a shared per-cell MLP                (Zaheer et al. 2017)
  gatedmil  : gated attention multiple-instance learning           (Ilse et al. 2018)
  pma       : pooling by multi-head attention with one seed vector (Lee et al. 2019)

No autograd framework is available in this environment, so every gradient is
derived analytically and checked against central finite differences
(`python poolnet.py selftest`).
"""
import numpy as np

def _seg_mean(A, starts, counts):
    return np.add.reduceat(A, starts, axis=0) / counts[:, None]

def sigmoid(x):
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    e = np.exp(x[~pos])
    out[~pos] = e / (1.0 + e)
    return out

def seg_softmax(e, starts, counts):
    e = np.atleast_2d(e.T).T if e.ndim == 1 else e
    mx = np.maximum.reduceat(e, starts, axis=0)
    ex = np.exp(e - np.repeat(mx, counts, axis=0))
    den = np.add.reduceat(ex, starts, axis=0)
    return ex / np.repeat(den, counts, axis=0)

def init_params(kind, K, H, rng, heads=4):
    s = 1.0 / np.sqrt(K)
    p = {"W1": rng.normal(0, s, (K, H)), "b1": np.zeros(H)}
    if kind == "deepsets":
        p["w2"] = np.zeros(H); p["b2"] = np.zeros(1)
    elif kind == "gatedmil":
        D = max(8, H // 2)
        t = 1.0 / np.sqrt(H)
        p["V"] = rng.normal(0, t, (H, D)); p["U"] = rng.normal(0, t, (H, D))
        p["wa"] = rng.normal(0, 1.0 / np.sqrt(D), D)
        p["w2"] = np.zeros(H); p["b2"] = np.zeros(1)
    elif kind == "pma":
        assert H % heads == 0
        p["q"] = rng.normal(0, 1.0 / np.sqrt(H // heads), (heads, H // heads))
        p["w2"] = np.zeros(H); p["b2"] = np.zeros(1)
    else:
        raise ValueError(kind)
    return p

def forward(kind, p, X, starts, counts, heads=4, cache=False):
    Z1 = X @ p["W1"] + p["b1"]
    A1 = np.maximum(Z1, 0.0)
    if kind == "deepsets":
        pooled = _seg_mean(A1, starts, counts)
        extra = {}
    elif kind == "gatedmil":
        T = np.tanh(A1 @ p["V"]); S = sigmoid(A1 @ p["U"])
        M = T * S
        a = seg_softmax(M @ p["wa"], starts, counts)[:, 0]
        pooled = np.add.reduceat(A1 * a[:, None], starts, axis=0)
        extra = {"T": T, "S": S, "M": M, "a": a}
    else:
        Hh = p["q"].shape[1]
        Ah = A1.reshape(len(A1), heads, Hh)
        e = np.einsum("nhd,hd->nh", Ah, p["q"]) / np.sqrt(Hh)
        a = seg_softmax(e, starts, counts)
        pooled = np.add.reduceat(Ah * a[:, :, None], starts, axis=0).reshape(len(starts), -1)
        extra = {"a": a, "Ah": Ah, "Hh": Hh}
    logit = pooled @ p["w2"] + p["b2"][0]
    if cache:
        return logit, {"Z1": Z1, "A1": A1, "pooled": pooled, **extra}
    return logit

def backward(kind, p, X, starts, counts, cache, glogit, heads=4):
    A1, pooled = cache["A1"], cache["pooled"]
    g = {k: np.zeros_like(v) for k, v in p.items()}
    g["w2"] = pooled.T @ glogit
    g["b2"] = np.array([glogit.sum()])
    dpooled = np.outer(glogit, p["w2"])
    if kind == "deepsets":
        dA1 = np.repeat(dpooled / counts[:, None], counts, axis=0)
    elif kind == "gatedmil":
        T, S, M, a = cache["T"], cache["S"], cache["M"], cache["a"]
        dP = np.repeat(dpooled, counts, axis=0)
        dA1 = dP * a[:, None]
        da = np.einsum("nh,nh->n", dP, A1)
        sda = np.repeat(np.add.reduceat(a * da, starts), counts)
        de = a * (da - sda)
        dM = np.outer(de, p["wa"])
        g["wa"] = M.T @ de
        dT = dM * S; dS = dM * T
        dZt = dT * (1.0 - T ** 2)
        dZs = dS * S * (1.0 - S)
        g["V"] = A1.T @ dZt
        g["U"] = A1.T @ dZs
        dA1 = dA1 + dZt @ p["V"].T + dZs @ p["U"].T
    else:
        a, Ah, Hh = cache["a"], cache["Ah"], cache["Hh"]
        dP = np.repeat(dpooled.reshape(len(starts), heads, Hh), counts, axis=0)
        dAh = dP * a[:, :, None]
        da = np.einsum("nhd,nhd->nh", dP, Ah)
        sda = np.repeat(np.add.reduceat(a * da, starts, axis=0), counts, axis=0)
        de = a * (da - sda) / np.sqrt(Hh)
        g["q"] = np.einsum("nh,nhd->hd", de, Ah)
        dAh = dAh + de[:, :, None] * p["q"][None, :, :]
        dA1 = dAh.reshape(len(A1), -1)
    dZ1 = dA1 * (cache["Z1"] > 0)
    g["W1"] = X.T @ dZ1
    g["b1"] = dZ1.sum(axis=0)
    return g

def bce_grad(logit, y, w):
    p = sigmoid(logit)
    return w * (p - y), float(np.sum(w * (np.logaddexp(0, logit) - y * logit)))

def train(kind, X, starts, counts, y, K, H, epochs, lr, wd, seed, heads=4, verbose=False):
    rng = np.random.default_rng(seed)
    p = init_params(kind, K, H, rng, heads)
    w = np.where(y == 1, 0.5 / max(y.mean(), 1e-9), 0.5 / max(1 - y.mean(), 1e-9))
    m = {k: np.zeros_like(v) for k, v in p.items()}
    v = {k: np.zeros_like(v) for k, v in p.items()}
    b1, b2, eps = 0.9, 0.999, 1e-8
    for t in range(1, epochs + 1):
        logit, cache = forward(kind, p, X, starts, counts, heads, cache=True)
        gl, loss = bce_grad(logit, y, w)
        g = backward(kind, p, X, starts, counts, cache, gl, heads)
        for k in p:
            if k not in ("b1", "b2"):
                g[k] = g[k] + wd * p[k]
            m[k] = b1 * m[k] + (1 - b1) * g[k]
            v[k] = b2 * v[k] + (1 - b2) * g[k] ** 2
            p[k] = p[k] - lr * (m[k] / (1 - b1 ** t)) / (np.sqrt(v[k] / (1 - b2 ** t)) + eps)
        if verbose and t % max(1, epochs // 5) == 0:
            print("  epoch %d loss %.4f" % (t, loss / len(y)))
    return p

def predict(kind, p, X, starts, counts, heads=4):
    return sigmoid(forward(kind, p, X, starts, counts, heads))

def _selftest():
    rng = np.random.default_rng(0)
    K, H, heads = 6, 8, 2
    counts = np.array([4, 3, 5])
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    X = rng.normal(size=(counts.sum(), K))
    y = np.array([1.0, 0.0, 1.0])
    w = np.ones(3)
    ok = True
    for kind in ("deepsets", "gatedmil", "pma"):
        p = init_params(kind, K, H, np.random.default_rng(1), heads)
        for k in p:
            p[k] = p[k] + rng.normal(0, 0.3, p[k].shape)
        logit, cache = forward(kind, p, X, starts, counts, heads, cache=True)
        gl, _ = bce_grad(logit, y, w)
        g = backward(kind, p, X, starts, counts, cache, gl, heads)
        worst = 0.0
        for k in p:
            flat = p[k].ravel()
            for i in range(min(len(flat), 12)):
                h = 1e-6
                old = flat[i]
                flat[i] = old + h
                _, lp = bce_grad(forward(kind, p, X, starts, counts, heads), y, w)
                flat[i] = old - h
                _, lm = bce_grad(forward(kind, p, X, starts, counts, heads), y, w)
                flat[i] = old
                num = (lp - lm) / (2 * h)
                ana = g[k].ravel()[i]
                worst = max(worst, abs(num - ana) / max(1.0, abs(num)))
        print("%-10s max relative gradient error = %.2e %s" % (kind, worst, "OK" if worst < 1e-5 else "FAIL"))
        ok = ok and worst < 1e-5
    print("SELFTEST", "PASS" if ok else "FAIL")
    return ok

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        raise SystemExit(0 if _selftest() else 1)
