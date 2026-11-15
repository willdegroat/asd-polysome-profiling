"""Check individual panel dimensions and approved pixel hashes."""

import argparse
import hashlib
import json
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    """Verify rendered panels against the approved source images."""
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--output', type = Path, default = ROOT / 'figures/generated')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'manuscript/panels.json').read_text())
    for record in manifest:
        path = args.output / record['panel'] / record['file']
        with Image.open(path) as image:
            assert list(image.size) == record['pixels'], path
            digest = hashlib.sha256(image.convert('RGB').tobytes()).hexdigest()
        assert digest == record['rgbSHA256'], (
            f'Pixels differ from the approved rendering: {path}. '
            'Check plotting versions and the Arial font installation.'
        )
        assert path.with_suffix('.pdf').exists(), path
        print(f'{record["panel"]}: dimensions and approved pixels match.')


if __name__ == '__main__':
    main()
