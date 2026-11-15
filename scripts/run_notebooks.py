import argparse
import os
import sys
from pathlib import Path
from time import perf_counter

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_DIR = ROOT / 'notebooks'
NOTEBOOKS = [path.name for path in sorted(NOTEBOOK_DIR.glob('*.ipynb'))]

parser = argparse.ArgumentParser()
parser.add_argument('notebooks', nargs = '*', choices = NOTEBOOKS)
args = parser.parse_args()
selected = args.notebooks or NOTEBOOKS

# The kernelspec uses `python`; prepend this environment to resolve it exactly.
python_dir = str(Path(sys.executable).resolve().parent)
os.environ['PATH'] = python_dir + os.pathsep + os.environ['PATH']

for name in selected:
    path = NOTEBOOK_DIR / name
    notebook = nbformat.read(path, as_version = 4)
    start = perf_counter()
    client = NotebookClient(
        notebook,
        timeout = None,
        kernel_name = 'python3',
        resources = {'metadata': {'path': str(ROOT)}},
        allow_errors = False,
    )
    client.execute()
    nbformat.write(notebook, path)
    print(f'{name}: {perf_counter() - start:.1f} seconds', flush = True)
