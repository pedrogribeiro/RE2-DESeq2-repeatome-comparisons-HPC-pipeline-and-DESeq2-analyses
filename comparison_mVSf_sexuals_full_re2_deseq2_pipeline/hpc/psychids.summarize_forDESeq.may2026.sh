#!/bin/bash
#PBS -N psychids.summarize_forDESeq.mVSf_sexuals
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=200gb
#PBS -l walltime=00:59:00
#PBS -o /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/prepareForDESeq2/mVSf_sexuals/out/summarize_forDESeq.txt
#PBS -e /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/prepareForDESeq2/mVSf_sexuals/err/summarize_forDESeq.txt

# Collapses per-base coverage (step6 output) into a cluster x sample count matrix for DESeq2.
# Cluster ID = contig name up to the first "_" (CLxxx_ prefix written when the mapping reference was built).

##### per-comparison knob #####
Comparison="mVSf_sexuals"

##### live log (plzen1) #####
LogFile="/storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/prepareForDESeq2/${Comparison}/liveLog/summarize_forDESeq.log"
mkdir -p "$(dirname "${LogFile}")"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LogFile}"; }

##### paths (brno12) #####
baseDir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/mapping/${Comparison}/out"            # INPUT: PSY_*/<Base>_PerBaseCov_sorted.csv
outDir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/prepareForDESeq2/${Comparison}/out"    # OUTPUT
outputMatrix="${outDir}/${Comparison}_cluster_matrix.tsv"
tmpDir="${outDir}/tmp_cluster_counts"
mkdir -p "${tmpDir}"

log "Job started"
log "Input directory: ${baseDir}"
log "Output matrix: ${outputMatrix}"

##### scratch setup #####
trap 'clean_scratch' TERM EXIT
cd "${SCRATCHDIR}" || exit 1
export TMPDIR="${SCRATCHDIR}"

##### step 1 - per-sample cluster counts #####
log "Step 1: generating per-sample cluster counts"

TotalSamples=$(ls -d "${baseDir}"/PSY_* 2>/dev/null | wc -l)
log "Sample directories found: ${TotalSamples}"

Counter=0
for sampleDir in "${baseDir}"/PSY_*; do
    Base=$(basename "${sampleDir}")
    input="${sampleDir}/${Base}_PerBaseCov_sorted.csv"
    localInput="${SCRATCHDIR}/${Base}_PerBaseCov_sorted.csv"
    output="${tmpDir}/${Base}.tsv"

    Counter=$((Counter + 1))
    log "Processing sample ${Counter}/${TotalSamples}: ${Base}"

    [ ! -f "${input}" ] && log "WARNING: input not found for ${Base}, skipping." && continue

    cp "${input}" "${localInput}"

    # sum column 3 (per-base depth) per cluster; cluster = column 1 truncated at first "_"
    awk '
    {
        cluster = $1
        sub(/_.*/, "", cluster)
        counts[cluster] += $3
    }
    END {
        for (c in counts) {
            print c "\t" counts[c]
        }
    }
    ' "${localInput}" | sort -V > "${output}"

    rm "${localInput}"    # free scratch per sample (otherwise every CSV accumulates in scratch)

    log "Done: ${Base}"
done

log "Step 1 complete"

##### step 2 - build matrix #####
log "Step 2: building cluster count matrix"

# full cluster list across all samples
cat "${tmpDir}"/*.tsv | cut -f1 | sort -V | uniq > "${tmpDir}/clusters.txt"

TotalClusters=$(wc -l < "${tmpDir}/clusters.txt")
log "Total clusters identified: ${TotalClusters}"

# header
echo -ne "Cluster" > "${outputMatrix}"
for f in "${tmpDir}"/*.tsv; do
    sample=$(basename "${f}" .tsv)
    echo -ne "\t${sample}" >> "${outputMatrix}"
done
echo >> "${outputMatrix}"

log "Matrix header written"

# fill matrix (missing cluster in a sample = 0)
while read -r cluster; do
    echo -ne "${cluster}" >> "${outputMatrix}"
    for f in "${tmpDir}"/*.tsv; do
        val=$(grep -w "^${cluster}" "${f}" | cut -f2)
        if [ -z "${val}" ]; then
            val=0
        fi
        echo -ne "\t${val}" >> "${outputMatrix}"
    done
    echo >> "${outputMatrix}"
done < "${tmpDir}/clusters.txt"

log "Step 2 complete. Matrix written to ${outputMatrix}"
log "Job finished"
