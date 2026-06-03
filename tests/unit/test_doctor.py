"""Unit tests for ssign_app/scripts/doctor.py — Resources section.

Other doctor sections (Python packages, External binaries, Databases,
Model weights) are exercised by tests/integration/test_doctor_*.py
against a real install. Here we only pin the new informational
Resources section so HPC users can spot scheduler-throttle bugs at
a glance.
"""

import io


class TestReportResources:
    def test_emits_cpu_ram_gpu_and_plme_lines(self, monkeypatch):
        from ssign_app.scripts import doctor

        # Pin the underlying detectors so the test is host-independent.
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_cpu_count", lambda: 8)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_ram_gb", lambda: 32.0)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.auto_batch_size_from_vram", lambda: 16)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.probe_cuda_device", lambda: ("FAKE-GPU", 24.0))
        monkeypatch.setattr("os.cpu_count", lambda: 16)

        buf = io.StringIO()
        doctor.report_resources(buf)
        out = buf.getvalue()

        assert "Resources" in out
        assert "effective: 8" in out and "host: 16" in out
        assert "32.0 GB" in out
        assert "FAKE-GPU, 24.0 GiB VRAM" in out
        assert "auto batch size: 16" in out

    def test_flags_scheduler_throttle(self, monkeypatch):
        # When effective_cpu_count < host cpu_count, doctor should call
        # out the gap — that's the bug we most want HPC users to notice.
        from ssign_app.scripts import doctor

        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_cpu_count", lambda: 4)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_ram_gb", lambda: 16.0)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.auto_batch_size_from_vram", lambda: 4)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.probe_cuda_device", lambda: (None, None))
        monkeypatch.setattr("os.cpu_count", lambda: 32)

        buf = io.StringIO()
        doctor.report_resources(buf)
        assert "scheduler is restricting" in buf.getvalue()

    def test_omits_throttle_note_when_full_allocation(self, monkeypatch):
        from ssign_app.scripts import doctor

        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_cpu_count", lambda: 16)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_ram_gb", lambda: 32.0)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.auto_batch_size_from_vram", lambda: 4)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.probe_cuda_device", lambda: (None, None))
        monkeypatch.setattr("os.cpu_count", lambda: 16)

        buf = io.StringIO()
        doctor.report_resources(buf)
        assert "scheduler is restricting" not in buf.getvalue()
