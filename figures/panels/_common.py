"""Locate publication inputs and individual panel output directories."""

import os
from pathlib import Path

from matplotlib import font_manager

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DATA = PROJECT_ROOT / 'manuscript/source_data'
ASSETS = PROJECT_ROOT / 'manuscript/assets'


def panel_output(label: str) -> Path:
    """Create one panel directory and require the approved typeface."""
    try:
        font_manager.findfont('Arial', fallback_to_default = False)
    except ValueError as error:
        raise RuntimeError(
            'Arial is required to reproduce the final panels. '
            'Install Arial before rendering.'
        ) from error
    root = Path(
        os.environ.get(
            'MANUSCRIPT_FIGURE_OUTPUT',
            PROJECT_ROOT / 'figures/generated',
        )
    )
    output = root / label
    output.mkdir(parents = True, exist_ok = True)
    return output
