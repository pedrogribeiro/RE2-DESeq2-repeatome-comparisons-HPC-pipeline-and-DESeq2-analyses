#!/bin/bash
# concatenateForClustering.sh - plain bash, run on the login node from this folder (no PBS needed).
# Concatenates the per-sample interleaved FASTAs from readProcessing into ONE file for RE2 clustering.
# Sample identity survives the concatenation: readProcessing stamps every read name as @PSY_xxxx/1 or /2,
# and RE2 recovers the sample from that prefix (--prefix_length 6).

set -euo pipefail

##### per-comparison knob #####
Comparison="geo"

##### paths (brno12) #####
Outdir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/readProcessing/${Comparison}/out"
Concat="${Outdir}/PSY${Comparison}Concatenated.fasta"

# never overwrite an existing concatenated FASTA (the real clustering input)
if [ -e "${Concat}" ]; then
    echo "ERROR: ${Concat} already exists, not overwriting." >&2
    exit 1
fi

# concatenate in samplesList.txt order; set -e aborts if any per-sample FASTA is missing
while IFS= read -r Base; do
    cat "${Outdir}/${Base}_12.fasta" >> "${Concat}"
done < samplesList.txt

echo "Written: ${Concat}"
