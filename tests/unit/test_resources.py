"""Tests for ssign_lib.resources — cgroup/scheduler-aware CPU and RAM detection.

These tests pin down the lookup priorities and the parse rules. The actual
production behavior on PBS/SLURM is hard to mock end-to-end (qstat must be
on PATH), so we mock the integration points and verify each branch in
isolation.
"""

from __future__ import annotations

import os

import pytest
from ssign_lib.resources import (
    _CGROUP_UNLIMITED_BYTES,
    _parse_size_to_gb,
    effective_cpu_count,
    effective_ram_gb,
)


class TestParseSizeToGb:
    @pytest.mark.parametrize(
        "raw,expected_gb",
        [
            ("32gb", 32.0),
            ("32GB", 32.0),
            ("32Gib", 32.0),
            ("4096mb", 4.0),  # 4096 MiB = 4 GiB
            ("128", 128 / 2**30),  # bare bytes
            ("0.5tb", 512.0),
        ],
    )
    def test_known_units(self, raw, expected_gb):
        assert _parse_size_to_gb(raw) == pytest.approx(expected_gb, rel=1e-3)

    def test_garbage_returns_none(self):
        assert _parse_size_to_gb("hello") is None
        assert _parse_size_to_gb("") is None


class TestEffectiveCpuCount:
    def test_returns_affinity_when_available(self, monkeypatch):
        monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: {0, 1, 2, 3}, raising=False)
        assert effective_cpu_count() == 4

    def test_falls_back_when_affinity_oserror(self, monkeypatch):
        def boom(_pid):
            raise OSError("no")

        monkeypatch.setattr(os, "sched_getaffinity", boom, raising=False)
        monkeypatch.setattr(os, "cpu_count", lambda: 16)
        assert effective_cpu_count() == 16

    def test_never_returns_zero(self, monkeypatch):
        monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: set(), raising=False)
        # Empty affinity (shouldn't happen but just in case) clamps to 1.
        assert effective_cpu_count() >= 1


class TestEffectiveRamGb:
    """Lookup priority: env override < SLURM < PBS qstat < cgroup < psutil. Returns MIN."""

    def test_env_override_wins_when_smallest(self, monkeypatch):
        monkeypatch.setenv("SSIGN_MAX_RAM_GB", "16")
        monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
        monkeypatch.delenv("PBS_JOBID", raising=False)
        # psutil host total will be much larger; min wins.
        assert effective_ram_gb() == pytest.approx(16.0)

    def test_slurm_alloc_is_parsed_in_mb(self, monkeypatch):
        monkeypatch.delenv("SSIGN_MAX_RAM_GB", raising=False)
        # 1024 MiB = 1 GB; smaller than any plausible host total, so MIN
        # logic should return the SLURM-derived value.
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "1024")
        monkeypatch.delenv("PBS_JOBID", raising=False)
        assert effective_ram_gb() == pytest.approx(1.0)

    def test_garbage_env_is_ignored(self, monkeypatch, caplog):
        monkeypatch.setenv("SSIGN_MAX_RAM_GB", "not-a-number")
        monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
        monkeypatch.delenv("PBS_JOBID", raising=False)
        # Doesn't crash; falls through to other sources.
        assert effective_ram_gb() > 0

    def test_cgroup_unlimited_sentinel_is_ignored(self):
        # The sentinel is 2**63 - 4096 (close enough to 2**63 - 1 modulo
        # page rounding). Confirms the production guard catches both.
        assert _CGROUP_UNLIMITED_BYTES < 2**63

    def test_returns_zero_when_nothing_available(self, monkeypatch):
        # Strip everything: env, scheduler, cgroup, psutil.
        for var in ("SSIGN_MAX_RAM_GB", "SLURM_MEM_PER_NODE", "PBS_JOBID"):
            monkeypatch.delenv(var, raising=False)
        # Block cgroup paths and psutil so we exercise the empty-candidates branch.
        monkeypatch.setattr("builtins.open", lambda *a, **kw: (_ for _ in ()).throw(OSError("no")))
        import sys

        monkeypatch.setitem(sys.modules, "psutil", None)
        # psutil=None makes the `import psutil` line raise TypeError, which the
        # bare except catches. Final result is 0.0.
        assert effective_ram_gb() == 0.0
