"""Shared failure-log helper for tool wrappers.

External-tool wrappers (run_bakta, run_eggnog, run_interproscan, ...) all
need to surface diagnostic info when their subprocess exits non-zero.
The naive pattern -- ``logger.error(result.stderr[:500])`` -- drops two
sources of evidence:

1. Some tools write their real error to stdout (Java stack traces from
   InterProScan), not stderr. Logging only stderr loses it.
2. CX3 compute-node ``$TMPDIR`` is auto-cleaned on job exit, so anything
   written there evaporates before the user can read it.

``dump_failure_log`` writes a single sidecar log file in the caller's
output directory, captures both streams in full, and returns a
``RuntimeError`` pre-populated with the path so the message itself is
the post-mortem entry point.

Usage::

    if result.returncode != 0:
        raise dump_failure_log("Bakta", result, cmd, output_dir)
"""

import logging
import os
import subprocess
from typing import Sequence

logger = logging.getLogger(__name__)

_LOG_STREAM_PREVIEW_CHARS = 500


def dump_failure_log(
    tool_name: str,
    result: subprocess.CompletedProcess,
    cmd: Sequence[str],
    output_dir: str,
) -> RuntimeError:
    """Persist a failed subprocess's full output and return a RuntimeError.

    The caller is expected to ``raise`` the returned exception. We don't
    raise inside the helper so call sites read as ``raise dump_failure_log(...)``
    rather than wrapping the side-effect in a separate function.

    Args:
        tool_name: Human-readable label for log messages and filename
            (e.g. ``"Bakta"``, ``"InterProScan"``). Spaces / dashes
            become underscores when used in the filename.
        result: CompletedProcess returned by ``subprocess.run(...,
            capture_output=True, text=True)``.
        cmd: The command list, recorded in the log for reproducibility.
        output_dir: Directory to write the sidecar log into. Must exist.

    Returns:
        ``RuntimeError`` with the exit code and absolute log path in its
        message, ready to ``raise``.
    """
    safe_name = tool_name.lower().replace(" ", "_").replace("-", "_")
    log_path = os.path.join(output_dir, f"{safe_name}_failure.log")
    try:
        with open(log_path, "w", encoding="utf-8") as f:
            f.write(f"{tool_name} exit code: {result.returncode}\n")
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"\n----- STDOUT -----\n{result.stdout or '(empty)'}\n")
            f.write(f"\n----- STDERR -----\n{result.stderr or '(empty)'}\n")
    except OSError as e:
        # Don't let a write failure mask the original tool failure --
        # the user needs to know the tool died, even if we couldn't
        # persist the details. Output_dir on RDS can hit quota or
        # permission errors mid-run.
        logger.warning("Could not write %s failure log to %s: %s", tool_name, log_path, e)
    logger.error(
        "%s exited with code %s. Full streams in %s\n  stdout (first %d): %s\n  stderr (first %d): %s",
        tool_name,
        result.returncode,
        log_path,
        _LOG_STREAM_PREVIEW_CHARS,
        (result.stdout or "(empty)")[:_LOG_STREAM_PREVIEW_CHARS],
        _LOG_STREAM_PREVIEW_CHARS,
        (result.stderr or "(empty)")[:_LOG_STREAM_PREVIEW_CHARS],
    )
    return RuntimeError(f"{tool_name} exit code {result.returncode}; full output in {log_path}")
