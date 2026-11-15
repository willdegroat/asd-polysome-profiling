"""Check publication archive integrity and restoration boundaries."""

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/restore_results.py'
SPEC = importlib.util.spec_from_file_location('restore_results', SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def make_archive(root: Path, *, destination: str = 'results/example.tsv') -> Path:
    """Create one deterministic archive for restoration contract checks."""
    data = b'geneId\tvalue\nexample\t1\n'
    packed = gzip.compress(data, mtime = 0)
    archive = root / 'published/example.tsv.gz'
    archive.parent.mkdir(parents = True)
    archive.write_bytes(packed)
    record = {
        'path': destination,
        'archive': 'published/example.tsv.gz',
        'bytes': len(data),
        'sha256': hashlib.sha256(data).hexdigest(),
        'archiveSHA256': hashlib.sha256(packed).hexdigest(),
    }
    (archive.parent / 'manifest.json').write_text(json.dumps({'files': [record]}))
    return archive


def test_restore_is_idempotent(tmp_path):
    make_archive(tmp_path)
    assert MODULE.restore(tmp_path) == 1
    assert MODULE.restore(tmp_path) == 0


def test_restore_rejects_modified_results(tmp_path):
    make_archive(tmp_path)
    MODULE.restore(tmp_path)
    target = tmp_path / 'results/example.tsv'
    target.write_text('new analysis\n')
    with pytest.raises(ValueError, match = 'differs from the frozen copy'):
        MODULE.restore(tmp_path)
    assert target.read_text() == 'new analysis\n'


def test_explicit_overwrite_restores_original_bytes(tmp_path):
    make_archive(tmp_path)
    MODULE.restore(tmp_path)
    target = tmp_path / 'results/example.tsv'
    target.write_text('new analysis\n')
    assert MODULE.restore(tmp_path, overwrite = True) == 1
    assert target.read_bytes() == b'geneId\tvalue\nexample\t1\n'


def test_corrupted_archive_is_rejected(tmp_path):
    archive = make_archive(tmp_path)
    archive.write_bytes(archive.read_bytes() + b'changed')
    with pytest.raises(ValueError, match = 'Archive checksum mismatch'):
        MODULE.restore(tmp_path)
    assert not (tmp_path / 'results/example.tsv').exists()


def test_paths_cannot_escape_checkout(tmp_path):
    make_archive(tmp_path, destination = '../outside.tsv')
    with pytest.raises(ValueError, match = 'inside the repository'):
        MODULE.restore(tmp_path)
