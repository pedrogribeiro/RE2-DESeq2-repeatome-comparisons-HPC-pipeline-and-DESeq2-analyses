#!/bin/bash

# Generate one PBS script per sample for mapping processed reads
# against the RepeatExplorer-derived mapping reference.
# Sample names are read from a list after: done < samplesList.txt

# input/output locations
ReadsDir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/readProcessing/mVSf_sexuals/out/intermediateOutputs/filtered"
IndexDir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/mapping/geo/out/refBuild"
ReferenceBase="geoConcatenatedRef"
Outdir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/mapping/mVSf_sexuals/out"

# Create output root directory if non-existent
mkdir -p "${Outdir}"

# Iterate from samples list to create running script
while IFS= read -r Base; do

    # input processed reads
    Forward="${ReadsDir}/${Base}_filt_1.fastq.gz"
    Reverse="${ReadsDir}/${Base}_filt_2.fastq.gz"

    # sample-specific output directory
    SampleOutdir="${Outdir}/${Base}"
    mkdir -p "${SampleOutdir}"

    # Name of PBS script to be written
    Outfile="${Base}.geoReMapping.pbs"

    (
    cat << EndOfJobScript
#!/bin/bash
#PBS -N ${Base}.psychids.geoReMapping.apr2026
#PBS -l select=1:ncpus=2:mem=164gb:scratch_local=50gb
#PBS -l walltime=23:59:00
#PBS -o /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/mVSf_sexuals/out/${Base}.geoReMapping.txt
#PBS -e /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/mVSf_sexuals/err/${Base}.geoReMapping.txt

## set live-log ##

LOG="/storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/mVSf_sexuals/liveLog/${Base}.progress.log"

set -euo pipefail

check() {
    echo "\$(date '+%Y-%m-%d %H:%M:%S') OK: \$1" | tee -a "\${LOG}" >&2
}

trap 'clean_scratch' TERM EXIT

##### system specific definitions #####
cd \$SCRATCHDIR || exit 1
export TMPDIR=\$SCRATCHDIR

# Copy input files and prebuilt bowtie2 index to scratch
rsync -a "${Forward}" "\$SCRATCHDIR/"
rsync -a "${Reverse}" "\$SCRATCHDIR/"
rsync -a "${IndexDir}/" "\$SCRATCHDIR/"
check "rsync input files"

# Define local scratch copies explicitly
LocalForward="${Base}_filt_1.fastq.gz"
LocalReverse="${Base}_filt_2.fastq.gz"

# Map reads to reference
module add bowtie2
bowtie2 -N 1 -X 101 --threads 2 -x "${ReferenceBase}" -1 "\${LocalForward}" -2 "\${LocalReverse}" -S "${Base}.sam"
module unload bowtie2
check "bowtie2 mapping"

# Samtools processing
module add samtools
samtools view -b -o "${Base}.bam" "${Base}.sam"
samtools view -F 256 -b "${Base}.bam" > "${Base}_primary.bam"
samtools view -H "${Base}_primary.bam" | grep @SQ | cut -d: -f2-3 | sed -e 's/\sLN:/\t/g' > "GenomeFile.csv"
samtools sort -m 3G -@ 2 -o "${Base}_prim_sorted.bam" "${Base}_primary.bam"
module unload samtools
check "samtools processing"

# Coverage
module add bedtools2
bedtools genomecov -d -ibam "${Base}_prim_sorted.bam" -g "GenomeFile.csv" > "${Base}_PerBaseCov_sorted.csv"
module unload bedtools2
check "bedtools coverage"

# Remove bulky intermediate files
rm "${Base}.bam"
rm "${Base}_primary.bam"
rm "${ReferenceBase}.1.bt2"
rm "${ReferenceBase}.2.bt2"
rm "${ReferenceBase}.3.bt2"
rm "${ReferenceBase}.4.bt2"
rm "${ReferenceBase}.fasta"
rm "${ReferenceBase}.rev.1.bt2"
rm "${ReferenceBase}.rev.2.bt2"
rm "${Base}_filt_1.fastq.gz"
rm "${Base}_filt_2.fastq.gz"

check "remove intermediate files"

# Copy outputs back
if rsync -a "\$SCRATCHDIR/" "${SampleOutdir}/"; then
    echo "OK: copy back" >&2
else
    export CLEAN_SCRATCH=false
    echo "FAILED: copy back - outputs may be lost" >&2
    exit 1
fi
EndOfJobScript
    ) > "${Outfile}"

    chmod +x "${Outfile}"
    echo "Script written: ${Outfile}"

done < samplesList.txt

SubmissionScript="submit.geoReMapping.jobs.sh"

(
cat << 'EndOfSubmissionScript'
#!/bin/bash
EndOfSubmissionScript

while IFS= read -r Base; do
    echo "qsub ${Base}.geoReMapping.pbs"
done < samplesList.txt

) > "${SubmissionScript}"

chmod +x "${SubmissionScript}"
echo "Submission script written: ${SubmissionScript}"
