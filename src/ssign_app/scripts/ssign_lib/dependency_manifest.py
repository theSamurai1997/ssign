"""Single source of truth for ssign's runtime dependencies.

Read by two consumers:

- ``src/ssign_app/scripts/doctor.py`` — backs the ``ssign doctor`` CLI command
  that users run on their own machine to verify a fresh install.
- ``tests/integration/test_imports.py`` — CI test that fails the moment
  ``pyproject.toml`` drifts from what the code actually imports.

Every entry carries a ``tier`` so the consumer can filter by install tier
(``base`` is always required; ``extended`` adds the annotation tools; ``full``
adds the largest databases).

If you add a new tool or pip dep to ``pyproject.toml``, add the matching entry
here in the same commit. The integration test will fail loudly otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

Tier = Literal["base", "extended", "full"]


@dataclass(frozen=True)
class PythonDep:
    """A Python package that must be importable.

    ``module`` is what ``importlib.import_module`` is given; ``pip_name`` is
    the pip-install name (often different — ``pyyaml`` provides ``yaml``,
    ``mkl-service`` provides ``mkl``). ``symbols`` are dotted attribute paths
    on top of the imported module that downstream code reaches for via lazy
    ``from X import Y`` style — checking them catches "package present but
    the symbol we use has moved/been removed" upstream-API drift.
    """

    module: str
    pip_name: str
    symbols: tuple[str, ...] = ()
    tier: Tier = "base"
    note: str = ""


@dataclass(frozen=True)
class ExternalBinary:
    """An executable that must be discoverable on ``PATH`` (or via a configured-install env var).

    ``install_dir_env`` names an env var holding a directory where the binary
    might also live — e.g. InterProScan ships as a tarball that the user
    extracts somewhere, then sets ``SSIGN_INTERPROSCAN_PATH``. Doctor checks
    ``PATH`` first, then ``${install_dir_env}/<binary>``.
    """

    name: str
    binary: str
    install_hint: str
    tier: Tier = "extended"
    optional: bool = False
    install_dir_env: str = ""


@dataclass(frozen=True)
class DatabasePath:
    """An on-disk database directory.

    Doctor resolves the path in this order: ``$SSIGN_<NAME>`` env var, then
    the ``$TARGET`` recorded in ``~/.ssign/db_root`` by the install script,
    then the default ``~/.ssign/databases/...``. ``sentinel_file`` is one
    specific file inside the directory whose presence confirms the DB is
    set up, not just that an empty directory exists.
    """

    name: str
    env_var: str
    default_subpath: str
    sentinel_file: str
    install_hint: str
    tier: Tier = "extended"


@dataclass(frozen=True)
class ModelWeights:
    """A model checkpoint / weights bundle.

    ``under_db_root`` toggles between the two layouts ssign uses:
      - False (default): subpath is relative to ``~/.ssign`` — covers things
        that auto-download on first pipeline run (DeepSecE checkpoint).
      - True: subpath is relative to the fetch_databases.sh target — covers
        bundles fetched by the install script (PLM-Effector weights).
    """

    name: str
    default_subpath: str
    install_hint: str
    tier: Tier = "base"
    under_db_root: bool = False


# ---------------------------------------------------------------------------
# Python packages
# ---------------------------------------------------------------------------

PYTHON_DEPS: tuple[PythonDep, ...] = (
    # Base tier — installed by `pip install ssign`
    PythonDep("streamlit", "streamlit"),
    PythonDep("pandas", "pandas"),
    PythonDep("numpy", "numpy"),
    PythonDep("Bio", "biopython", symbols=("Bio.SeqIO", "Bio.SeqUtils.ProtParam.ProteinAnalysis")),
    PythonDep("matplotlib", "matplotlib"),
    PythonDep("seaborn", "seaborn"),
    PythonDep("scipy", "scipy"),
    PythonDep("requests", "requests"),
    PythonDep("yaml", "pyyaml"),
    PythonDep("jinja2", "jinja2"),
    PythonDep("networkx", "networkx"),
    PythonDep("taxopy", "taxopy"),
    PythonDep("pyhmmer", "pyhmmer"),
    PythonDep("pyrodigal", "pyrodigal"),
    PythonDep("macsylib", "macsyfinder", note="macsyfinder pip package ships the `macsylib` module"),
    PythonDep(
        "esm",
        "fair-esm",
        symbols=("esm.Alphabet", "esm.FastaBatchedDataset"),
    ),
    PythonDep("torch", "torch", symbols=("torch.cuda", "torch.serialization", "torch.utils.data.DataLoader")),
    PythonDep(
        "DeepSecE",
        "deepsece",
        symbols=("DeepSecE.model.EffectorTransformer",),
        note="lazy-imported from run_deepsece.py at runtime",
    ),
    # Extended tier — installed by `pip install ssign[extended]`
    PythonDep("bakta", "bakta", tier="extended"),
    PythonDep("goatools", "goatools", tier="extended"),
    PythonDep("obonet", "obonet", tier="extended"),
    PythonDep(
        "transformers",
        "transformers",
        symbols=(
            "transformers.T5Tokenizer",
            "transformers.T5EncoderModel",
            "transformers.AutoModel",
            "transformers.AutoTokenizer",
        ),
        tier="extended",
        note="used by PLM-Effector (T5*) + pLM-BLAST embedders (prottrans, hfautomodel)",
    ),
    PythonDep(
        "xgboost",
        "xgboost",
        symbols=("xgboost.XGBClassifier",),
        tier="extended",
        note="PLM-Effector ensemble.py",
    ),
    PythonDep(
        "google.protobuf",
        "protobuf",
        tier="extended",
        note="required by T5 SentencePiece tokenizer",
    ),
    PythonDep(
        "mkl",
        "mkl",
        tier="extended",
        note="`mkl` module is also provided by the `mkl-service` pip package; both are listed",
    ),
    PythonDep("h5py", "h5py", tier="extended", note="pLM-BLAST embedders/dataset.py"),
    PythonDep("tqdm", "tqdm", tier="extended"),
    PythonDep("fairscale", "fairscale", tier="extended", note="pLM-BLAST ESM FSDP wrapper"),
    PythonDep("numba", "numba", tier="extended", note="pLM-BLAST scripts/plmblast.py + alntools/numeric"),
    PythonDep("sentencepiece", "sentencepiece", tier="extended", note="T5Tokenizer backend"),
    PythonDep("accelerate", "accelerate", tier="extended", note="transformers HF model loaders"),
)


# ---------------------------------------------------------------------------
# External binaries (not pip-installable)
# ---------------------------------------------------------------------------

EXTERNAL_BINARIES: tuple[ExternalBinary, ...] = (
    ExternalBinary("Bakta", "bakta", "pip install ssign[extended]  # bakta ships as a pip package"),
    ExternalBinary(
        "EggNOG-mapper",
        "emapper.py",
        "separate conda env required (biopython conflict with bakta) — see docs/how-to/install.md § EggNOG",
    ),
    ExternalBinary("HH-suite — hhsearch", "hhsearch", "module load HH-suite, or conda install -c bioconda hhsuite"),
    ExternalBinary("HH-suite — hhblits", "hhblits", "module load HH-suite, or conda install -c bioconda hhsuite"),
    ExternalBinary("BLAST+", "blastp", "module load BLAST+, or apt install ncbi-blast+"),
    ExternalBinary(
        "InterProScan",
        "interproscan.sh",
        "see docs/how-to/install.md § InterProScan (Java + 30 GB tarball)",
        install_dir_env="SSIGN_INTERPROSCAN_PATH",
    ),
    ExternalBinary(
        "SignalP 6 (local)",
        "signalp6",
        "DTU portal — only required if --signalp-mode=local (default since 723a96f)",
        tier="base",
        optional=True,
    ),
    ExternalBinary(
        "DeepLocPro (local)",
        "deeplocpro",
        "DTU portal — only required if --deeplocpro-mode=local (default since 723a96f)",
        tier="base",
        optional=True,
    ),
)


# ---------------------------------------------------------------------------
# Database paths (default root: ~/.ssign/databases/ per scripts/fetch_databases.sh)
# ---------------------------------------------------------------------------
# default_subpath is relative to either $SSIGN_<NAME> (if set) or ~/.ssign/databases/

DATABASE_PATHS: tuple[DatabasePath, ...] = (
    DatabasePath(
        "Bakta DB",
        "SSIGN_BAKTA_DB",
        "bakta",
        "version.json",
        "bash scripts/fetch_databases.sh --tier base",
        tier="extended",
    ),
    DatabasePath(
        "EggNOG DB",
        "SSIGN_EGGNOG_DB",
        "eggnog",
        "eggnog.db",
        "bash scripts/fetch_databases.sh --tier extended",
        tier="extended",
    ),
    DatabasePath(
        "InterProScan DB",
        "SSIGN_INTERPROSCAN_PATH",
        "interproscan",
        "interproscan.properties",
        "bash scripts/fetch_databases.sh --tier extended",
        tier="extended",
    ),
    DatabasePath(
        "HH-suite Pfam",
        "SSIGN_HHSUITE_PFAM",
        "hhsuite/pfam",
        "pfam_a3m.ffdata",
        "bash scripts/fetch_databases.sh --tier extended",
        tier="extended",
    ),
    DatabasePath(
        "HH-suite PDB70",
        "SSIGN_HHSUITE_PDB70",
        "hhsuite/pdb70",
        "pdb70_a3m.ffdata",
        "bash scripts/fetch_databases.sh --tier extended",
        tier="extended",
    ),
    DatabasePath(
        "HH-suite UniRef30",
        "SSIGN_HHSUITE_UNICLUST",
        "hhsuite/uniref30",
        "UniRef30_2023_02_a3m.ffdata",
        "bash scripts/fetch_databases.sh --tier full",
        tier="full",
    ),
    DatabasePath(
        "pLM-BLAST ECOD70",
        "SSIGN_ECOD70_DB",
        "plm_blast",
        "ECOD70.csv",
        "bash scripts/fetch_databases.sh --tier extended",
        tier="extended",
    ),
    DatabasePath(
        "BLAST NR",
        "SSIGN_BLAST_NR",
        "blast_nr",
        "nr.pdb",
        "bash scripts/fetch_databases.sh --tier full",
        tier="full",
    ),
)


# ---------------------------------------------------------------------------
# Model weights (default root: ~/.ssign/models/, or ~/.ssign/databases/ for fetched bundles)
# ---------------------------------------------------------------------------

MODEL_WEIGHTS: tuple[ModelWeights, ...] = (
    ModelWeights(
        "DeepSecE checkpoint",
        "models/deepsece_checkpoint.pt",
        "auto-downloaded on first run, or `bash scripts/fetch_weights.sh`",
    ),
    ModelWeights(
        "PLM-Effector ensemble weights",
        "plm_effector_weights",
        "bash scripts/fetch_weights.sh (or scripts/fetch_databases.sh --tier extended)",
        tier="extended",
        under_db_root=True,
    ),
)


_TIER_ORDER: tuple[Tier, ...] = ("base", "extended", "full")


def _filter_by_tier(items, tier: Tier):
    upto = _TIER_ORDER[: _TIER_ORDER.index(tier) + 1]
    return tuple(i for i in items if i.tier in upto)


def deps_for_tier(tier: Tier) -> tuple[PythonDep, ...]:
    return _filter_by_tier(PYTHON_DEPS, tier)


def binaries_for_tier(tier: Tier) -> tuple[ExternalBinary, ...]:
    return _filter_by_tier(EXTERNAL_BINARIES, tier)


def databases_for_tier(tier: Tier) -> tuple[DatabasePath, ...]:
    return _filter_by_tier(DATABASE_PATHS, tier)


def weights_for_tier(tier: Tier) -> tuple[ModelWeights, ...]:
    return _filter_by_tier(MODEL_WEIGHTS, tier)
