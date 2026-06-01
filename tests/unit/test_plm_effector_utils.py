"""Tests for ssign_app.scripts.plm_effector.utils.resolve_autocast_dtype.

The helper maps user-facing dtype names ("bf16", "fp16", "fp32", aliases)
to the torch.dtype consumed by torch.autocast. None means "skip autocast",
which the call site uses to keep the model in its loaded precision.
"""

import sys

import pytest
import torch

# Make the script dir importable so `from plm_effector.utils import ...`
# resolves the same way it does in production (script subprocess on PATH).
_SCRIPT_DIR = (
    __import__("os")
    .path.dirname(__import__("os").path.abspath(__file__))
    .replace("/tests/unit", "/src/ssign_app/scripts")
)
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from plm_effector.utils import resolve_autocast_dtype  # noqa: E402


class TestResolveAutocastDtype:
    @pytest.mark.parametrize("name", ["fp32", "float32", "none", "", None])
    def test_fp32_aliases_return_none(self, name):
        assert resolve_autocast_dtype(name) is None

    @pytest.mark.parametrize("name", ["bf16", "bfloat16", "BF16", "  bfloat16  "])
    def test_bf16_aliases_return_bfloat16(self, name):
        assert resolve_autocast_dtype(name) is torch.bfloat16

    @pytest.mark.parametrize("name", ["fp16", "float16", "half", "FP16"])
    def test_fp16_aliases_return_float16(self, name):
        assert resolve_autocast_dtype(name) is torch.float16

    def test_torch_dtype_passes_through(self):
        assert resolve_autocast_dtype(torch.bfloat16) is torch.bfloat16
        assert resolve_autocast_dtype(torch.float16) is torch.float16

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError, match="Unknown PLM-Effector dtype"):
            resolve_autocast_dtype("int8")
