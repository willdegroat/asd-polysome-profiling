"""Render approved manuscript panels individually."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from PIL import Image
from reportlab.pdfgen import canvas

ROOT = Path(__file__).resolve().parents[1]
PANELS = [
    '1A',
    '1B',
    '1C',
    '1D',
    '1E',
    '1legend',
    '2A',
    '2B',
    '2C',
    '2D',
    '2E',
    '2F',
    '2G',
    '2H',
    '2I',
    '2legend',
    '3A',
    '3B',
    '3C',
    '3D',
]


def main() -> None:
    """Render selected panels without assembling full figures."""
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--panels', nargs = '+', choices = PANELS, default = PANELS)
    parser.add_argument(
        '--output', type = Path, default = ROOT / 'figures/generated'
    )
    args = parser.parse_args()
    if not (ROOT / 'results/polysomeMonosome.tsv').exists():
        raise FileNotFoundError('Run python scripts/restore_results.py first.')
    env = dict(os.environ)
    cache = args.output.resolve() / '.mplconfig'
    cache.mkdir(parents = True, exist_ok = True)
    env['MPLCONFIGDIR'] = str(cache)
    env['MANUSCRIPT_FIGURE_OUTPUT'] = str(args.output.resolve())
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for panel in args.panels:
        script = ROOT / 'figures/panels' / (panel.lower() + '.py')
        subprocess.run(
            [sys.executable, '-B', str(script)], cwd = ROOT, env = env, check = True
        )
        for image in (args.output / panel).glob('Figure_*.png'):
            with Image.open(image) as bitmap:
                width, height = (value * 72 / 600 for value in bitmap.size)
            pdf = canvas.Canvas(
                str(image.with_suffix('.pdf')), pagesize = (width, height)
            )
            pdf.drawImage(str(image), 0, 0, width = width, height = height)
            pdf.showPage()
            pdf.save()
        print(f'Rendered {panel} in {args.output / panel}', flush = True)


if __name__ == '__main__':
    main()
