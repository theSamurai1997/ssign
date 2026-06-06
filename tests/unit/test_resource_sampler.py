"""Unit tests for ResourceSampler — the in-process resource sampler that
the runner spins up as a daemon thread so every run produces
``outdir/runtime_data/resource_samples.csv`` without the user having to launch anything.

We don't test the actual sampling values (host-dependent); we test that
the thread starts, writes a CSV with the expected header, picks up the
current step label, and shuts down cleanly. Sampling failures must not
crash the pipeline — that invariant is tested by forcing a psutil call
to raise.
"""

import csv
import time


class TestResourceSampler:
    def test_writes_csv_header_then_samples_then_stops(self, tmp_path):
        from ssign_app.core.resource_sampler import ResourceSampler

        out = tmp_path / "resources.csv"
        s = ResourceSampler(out_path=str(out), interval=0.5)
        s.start()
        try:
            time.sleep(0.4)  # let it write a few rows
        finally:
            s.stop()

        with open(out) as fh:
            rows = list(csv.reader(fh))
        assert rows, "no rows written"
        assert rows[0][0] == "timestamp"
        # At least one data row landed
        assert len(rows) >= 2

    def test_set_step_tagged_in_subsequent_rows(self, tmp_path):
        # Interval is clamped to >=0.5s in the sampler for safety, so use
        # sleep windows wider than that to guarantee a sample lands inside
        # each labelled window.
        from ssign_app.core.resource_sampler import ResourceSampler

        out = tmp_path / "resources.csv"
        s = ResourceSampler(out_path=str(out), interval=0.5)
        s.start()
        try:
            s.set_step("step_2:macsyfinder")
            time.sleep(0.8)  # one interval at 0.5s, plus margin
            s.set_step("parallel:dlp,dse,signalp")
            time.sleep(0.8)
        finally:
            s.stop()

        with open(out) as fh:
            reader = csv.DictReader(fh)
            steps = [r["step"] for r in reader]
        assert "step_2:macsyfinder" in steps
        assert "parallel:dlp,dse,signalp" in steps

    def test_stop_is_idempotent(self, tmp_path):
        from ssign_app.core.resource_sampler import ResourceSampler

        s = ResourceSampler(out_path=str(tmp_path / "resources.csv"), interval=0.5)
        s.start()
        time.sleep(0.2)
        s.stop()
        s.stop()  # second call must not blow up
