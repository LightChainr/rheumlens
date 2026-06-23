from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.manifest, sep='\t', dtype=str).fillna('')
    statuses: list[str] = []
    observed_hashes: list[str] = []
    for row in table.itertuples(index=False):
        raw_path = str(row.path_or_source).strip()
        path = Path(raw_path).expanduser() if raw_path and not raw_path.startswith(('http://', 'https://')) else None
        if path is not None and path.exists() and path.is_file():
            statuses.append('ready')
            observed_hashes.append(sha256_file(path))
        elif path is not None and path.exists():
            statuses.append('ready_dir')
            observed_hashes.append('')
        else:
            statuses.append(str(row.status) if str(row.status) else 'missing')
            observed_hashes.append('')
    table['observed_status'] = statuses
    table['observed_sha256'] = observed_hashes
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    missing_required = table[(table['required'].str.lower() == 'yes') & ~table['observed_status'].isin(['ready', 'ready_dir'])]
    print(f'assets: {len(table)}')
    print(f'missing required: {len(missing_required)}')
    if len(missing_required):
        print(missing_required[['asset_id', 'category', 'cohort_or_model', 'path_or_source']].to_string(index=False))
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
