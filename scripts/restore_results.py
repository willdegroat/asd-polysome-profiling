"""Restore the frozen manuscript results without refitting models."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def restore(root: Path, *, overwrite: bool = False) -> int:
    """Restore checked archives and reject modified existing results."""
    root = root.resolve()
    manifest = json.loads((root / 'published/manifest.json').read_text())
    restored = 0
    for record in manifest['files']:
        source = (root / record['archive']).resolve()
        target = (root / record['path']).resolve()
        if not source.is_relative_to(root) or not target.is_relative_to(root):
            raise ValueError('Archive paths must stay inside the repository.')
        compressed = source.read_bytes()
        if hashlib.sha256(compressed).hexdigest() != record['archiveSHA256']:
            raise ValueError(f'Archive checksum mismatch: {record["archive"]}')
        if target.exists():
            current = hashlib.sha256(target.read_bytes()).hexdigest()
            if current == record['sha256']:
                continue
            if not overwrite:
                raise ValueError(
                    f'Result differs from the frozen copy: {record["path"]}. '
                    'Use a separate checkout, or --overwrite to restore it.'
                )
        data = gzip.decompress(compressed)
        if len(data) != record['bytes']:
            raise ValueError(f'Result size mismatch: {record["path"]}')
        if hashlib.sha256(data).hexdigest() != record['sha256']:
            raise ValueError(f'Result checksum mismatch: {record["path"]}')
        target.parent.mkdir(parents = True, exist_ok = True)
        target.write_bytes(data)
        restored += 1
    return restored


def main() -> None:
    """Restore the published archive into this checkout."""
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--overwrite', action = 'store_true')
    args = parser.parse_args()
    restored = restore(ROOT, overwrite = args.overwrite)
    print(f'Restored {restored} files; all archives passed SHA-256 checks.')


if __name__ == '__main__':
    main()
