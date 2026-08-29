"""Source-file fingerprinting: strong vs fast mode, sensitivity, map behaviour."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

from analysis.campaign import fingerprint as fp_mod
from analysis.campaign.fingerprint import fingerprint_file, fingerprint_map


def test_small_file_is_strong_full_sha256(tmp_path):
    p = tmp_path / "small.xtc"
    payload = b"hello trajectory" * 100
    p.write_bytes(payload)

    fp = fingerprint_file(p)
    assert fp.mode == "strong"
    assert fp.digest == hashlib.sha256(payload).hexdigest()
    assert fp.size_bytes == len(payload)
    assert fp.short == fp.digest[:12]


def test_strong_flag_forces_full_hash(tmp_path, monkeypatch):
    p = tmp_path / "big.xtc"
    payload = b"x" * 4096
    p.write_bytes(payload)
    monkeypatch.setattr(fp_mod, "STRONG_MAX_BYTES", 16)  # would otherwise be "fast"

    fp = fingerprint_file(p, strong=True)
    assert fp.mode == "strong"
    assert fp.digest == hashlib.sha256(payload).hexdigest()


def test_fast_mode_when_over_threshold(tmp_path, monkeypatch):
    monkeypatch.setattr(fp_mod, "STRONG_MAX_BYTES", 8)
    p = tmp_path / "traj.xtc"
    p.write_bytes(b"A" * 512)

    fp = fingerprint_file(p)
    assert fp.mode == "fast"
    assert fp.edge_bytes == fp_mod.FAST_EDGE_BYTES
    assert "not collision-proof" in fp.note
    # fast digest is NOT the plain content hash
    assert fp.digest != hashlib.sha256(b"A" * 512).hexdigest()


def test_fast_digest_changes_with_mtime(tmp_path, monkeypatch):
    monkeypatch.setattr(fp_mod, "STRONG_MAX_BYTES", 8)
    p = tmp_path / "traj.xtc"
    p.write_bytes(b"A" * 512)

    d1 = fingerprint_file(p).digest
    future = (p.stat().st_atime_ns + 10 ** 12, p.stat().st_mtime_ns + 10 ** 12)
    os.utime(p, ns=future)
    d2 = fingerprint_file(p).digest
    assert d1 != d2


def test_fast_digest_changes_with_size(tmp_path, monkeypatch):
    monkeypatch.setattr(fp_mod, "STRONG_MAX_BYTES", 8)
    p = tmp_path / "traj.xtc"
    p.write_bytes(b"A" * 512)
    d1 = fingerprint_file(p).digest
    p.write_bytes(b"A" * 1024)
    d2 = fingerprint_file(p).digest
    assert d1 != d2


def test_fast_digest_changes_with_edge_content(tmp_path, monkeypatch):
    monkeypatch.setattr(fp_mod, "STRONG_MAX_BYTES", 8)
    p = tmp_path / "traj.xtc"
    p.write_bytes(b"A" * 512)
    mtime = (p.stat().st_atime_ns, p.stat().st_mtime_ns)
    d1 = fingerprint_file(p).digest

    p.write_bytes(b"B" * 512)          # same size, different edge bytes
    os.utime(p, ns=mtime)             # pin mtime so only content differs
    d2 = fingerprint_file(p).digest
    assert d1 != d2


def test_fingerprint_map_skips_missing_and_none(tmp_path):
    real = tmp_path / "real.gro"
    real.write_text("data")
    out = fingerprint_map({
        "structure": real,
        "topology": tmp_path / "does_not_exist.top",
        "trajectory": None,
    })
    assert set(out) == {"structure"}
    assert out["structure"].mode == "strong"
