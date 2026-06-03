# Awk script consumed by scripts/cx3/submit_all_gpus.sh.
#
# Reads `pbsnodes -a` output on stdin and prints every distinct
# gpu_type advertised by a node whose state is usable (not down /
# offline / stale / state-unknown). Drops the literal "None" gpu_type
# (CX3 shorthand for CPU-only nodes) and any empty-value entry.
#
# Sort/dedup the output yourself: `awk -f this.awk | sort -u`.

function emit() {
    # PBS Pro node state is comma-separated like "free", "job-busy",
    # "job-busy,resv-exclusive", or "down,offline". Any "down"/
    # "offline"/"stale"/"state-unknown" marker means jobs sent there
    # won't start; treat the gpu_type as unusable.
    if (gpu != "" && gpu != "None" && state != "" \
        && state !~ /down|offline|stale|state-unknown/)
        print gpu
    gpu = ""
    state = ""
}

# Hostname line (no leading whitespace) marks a new node block.
/^[^ \t]/                                    { emit() }

/^[ \t]+state *=/                            { state = $0 }

/resources_available\.gpu_type *=/ {
    split($0, parts, /= */)
    v = parts[2]
    gsub(/^[ \t]+|[ \t]+$/, "", v)
    gpu = v
}

END { emit() }
