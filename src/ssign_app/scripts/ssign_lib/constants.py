"""Central constants and default thresholds for the ssign pipeline.

All configurable thresholds are defined here so that (a) the Python scripts
have sensible defaults when run standalone and (b) the Nextflow params
override them at invocation time via CLI arguments.
"""

# --- DeepLocPro ---
CONF_THRESHOLD = 0.8  # Minimum extracellular probability

# --- MacSyFinder ---
WHOLENESS_THRESHOLD = 0.8  # Minimum system completeness
REQUIRED_FRACTION_CORRECT = 0.8  # Fraction of SS components correctly localized

# --- Proximity ---
PROXIMITY_WINDOW = 3  # +/- N genes per SS component

# --- Structure quality ---
PLDDT_THRESHOLD = 70  # Minimum mean pLDDT for structure acceptance

# --- BLASTp ---
BLASTP_MIN_PIDENT = 80  # Minimum percent identity
BLASTP_MIN_QCOV = 80  # Minimum query coverage

# --- Tool subprocess timeout ---
# Generic upper bound for any external bioinformatics tool ssign shells
# out to (BLAST, IPS, EggNOG, Bakta, pLM-BLAST, DLP, SignalP). 4h is
# deliberately generous — pathological inputs (giant proteome, huge DB,
# slow disk) can stretch any of these tools well past the typical
# 10-60 min runtime. A timeout firing means the tool is wedged or the
# input is unrealistically large; give up rather than block the rest of
# the pipeline forever. Bump if you see legitimate timeouts.
# (HH-suite has its own per-stage timeouts a few sections below — those
# are tighter because they fire per-protein, not per-batch.)
TOOL_TIMEOUT_S = 14400  # 4h

# --- HH-suite ---
# HHR Prob (0-100) cutoff for keeping the top-1 hit per DB. Söding-lab
# guidance: ≥95 near-certain homolog, ≥50 worth considering. ssign default
# of 80 is permissive-but-meaningful. Filtering on Prob, not E-value.
HHSUITE_MIN_PROB = 80.0
# Three iterations balance sensitivity for remote homologs against
# runtime; ssign exists to surface distant homologs of secretion substrates.
HHBLITS_ITERATIONS = 3
# Per-protein subprocess timeouts. Empirical wall-time for a 1500-aa
# autotransporter against UniRef30 (~204 GB) with -n 3 iterations on a
# 4-thread CPU run is 30-60 min for hhblits; hhsearch against Pfam/PDB70
# is faster (~5-15 min) but capped at half the hhblits budget so total
# wall-time per protein stays bounded.
HHBLITS_TIMEOUT_S = 3600
HHSEARCH_TIMEOUT_S = 1800

# --- DTU API HTTP timeouts (run_deeplocpro.py + run_signalp.py) ---
# Three tiers matching what the DTU webserver actually does at each call.
# Submit is a multipart upload that the server queues — generous budget
# accounts for the upload itself plus the server's `?ajax=1` redirect dance.
# Status polls are a single JSON GET against the queue, so much shorter.
# Result download is a tar fetch, in between.
DTU_API_SUBMIT_TIMEOUT_S = 60
DTU_API_STATUS_TIMEOUT_S = 15
DTU_API_DOWNLOAD_TIMEOUT_S = 30

# --- map_gbff_to_bakta_cds.py ---
# Minimum coordinate-overlap fraction used when matching a Bakta CDS back
# to its GenBank counterpart. 0.8 catches near-identical predictions while
# rejecting spurious matches caused by gene-prediction tool disagreement.
MAP_GBFF_BAKTA_MIN_OVERLAP = 0.8

# --- T5aSS domain classification ---
MIN_PASSENGER_LENGTH = 100  # aa — below this = "minimal passenger"
LINKER_LENGTH = 30  # aa — alpha-helix linker between passenger and barrel

# --- System filtering ---
DEFAULT_EXCLUDED_SYSTEMS = ["Flagellum", "Tad", "T3SS"]

# --- T5SS per-component DLP rules ---
# Per-component biology rationale lives in cross_validate_predictions.py's
# module docstring (the consumer); this dict is the data.
#
# FRAGILE: gene_names are TXSScan v2.1 model IDs. If TXSScan renames a
# component (Pfam succession, model rebuild), this dict goes stale
# silently — affected subtypes revert to the strict ext-only rule and
# under-call passengers. If this breaks: pin TXSScan or update keys to
# match the new model names.
T5SS_COMPONENT_RULES = {
    ("T5aSS", "T5aSS_PF03797"): ("extracellular_prob", "outer_membrane_prob"),
    ("T5bSS", "T5bSS_translocator"): ("outer_membrane_prob",),
    ("T5cSS", "T5cSS_PF03895"): ("extracellular_prob", "outer_membrane_prob"),
}

# --- DeepSecE to MacSyFinder SS type mapping ---
# DeepSecE predicts broad types; MacSyFinder uses specific names.
DSE_TO_MACSYFINDER = {
    "T1SS": ["T1SS"],
    "T2SS": ["T2SS"],
    "T3SS": ["T3SS"],
    "T4SS": ["pT4SSt", "T4SS"],
    "T6SS": ["T6SSi", "T6SS"],
}

# --- PLM-Effector type to MacSyFinder SS type mapping ---
# PLM-Effector reports per-effector-type predictions (T1SE = T1-secreted-
# effector, etc.); the gate at proximity time needs to map those back to
# MacSyFinder SS-system names to verify the system actually exists in the
# genome (same cross-genome leakage class of bug DSE_TO_MACSYFINDER fixes).
PLME_TO_MACSYFINDER = {
    "T1SE": ["T1SS"],
    "T2SE": ["T2SS"],
    "T3SE": ["T3SS"],
    "T4SE": ["pT4SSt", "T4SS"],
    "T6SE": ["T6SSi", "T6SS"],
}


# --- Install-tier defaults ---
# Which optional tools default on at each install tier. The user's tier
# comes from --tier on the CLI, or ~/.ssign/tier (written by
# fetch_databases.sh at the end of a successful fetch), or falls back to
# "extended". Per-tool overrides (--skip-X / --no-skip-X) still win when
# the user is explicit.
#
# Mapping reflects "what does this tier ship that's also actually usable?":
#   - base: Bakta light DB + DeepSecE/PLM-Effector weights → those tools
#     default on. Annotation tools that need extended DBs (EggNOG ~50 GB,
#     IPS ~24 GB, ECOD70 ~24 GB) default off.
#   - extended: adds those three DBs → EggNOG/IPS/pLM-BLAST default on.
#     Pfam + PDB70 also get fetched, but HH-suite stays default-off (see
#     _EXTENDED_ADDS comment) until the wrapper degrades gracefully
#     without UniRef30.
#   - full: adds BLAST NR (390 GB) + HH-suite UniRef30 → BLASTp and
#     HH-suite default on at the full tier.
#
# `run_bakta` is governed by --use-input-annotations + Bakta-DB-present
# invariant separately; not in this table.
# Each tool listed exactly once at the tier it first becomes available.
# The full per-tier on/off map is built below.
_BASE_ENABLED = frozenset({"deeplocpro", "signalp", "deepsece", "plm_effector", "protparam"})
# HH-suite is NOT in extended because its hhblits MSA step needs UniRef30
# (~25 GB), which only ships with --tier full. The wrapper currently
# aborts when UniRef30 is missing rather than degrading to hhsearch-only;
# until that lands, default-on at extended would mislead users into
# thinking HH-suite is available when it can't actually run. Pfam +
# PDB70 still get fetched at extended for when the user wants to enable
# HH-suite manually (--no-skip-hhsuite) on a node with UniRef30 present.
_EXTENDED_ADDS = frozenset({"interproscan", "eggnog", "plmblast"})
_FULL_ADDS = frozenset({"blastp", "hhsuite"})

_TIER_ENABLED = {
    "base": _BASE_ENABLED,
    "extended": _BASE_ENABLED | _EXTENDED_ADDS,
    "full": _BASE_ENABLED | _EXTENDED_ADDS | _FULL_ADDS,
}
_ALL_TOOLS = _TIER_ENABLED["full"]

TIER_TOOL_DEFAULTS = {tier: {tool: (tool in enabled) for tool in _ALL_TOOLS} for tier, enabled in _TIER_ENABLED.items()}
DEFAULT_TIER = "extended"
