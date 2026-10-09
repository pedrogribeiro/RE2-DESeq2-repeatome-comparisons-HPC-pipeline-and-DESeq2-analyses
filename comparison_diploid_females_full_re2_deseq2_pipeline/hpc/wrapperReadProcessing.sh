#!/bin/bash

# Generate one PBS script per sample for RepeatExplorer read preprocessing
# extracting names from a list after done < samplesList.txt 
# This script: trims and filters reads, samples them so RepeatExplorer can deal with amount
# and finally converts them to interleaved FASTA. 

# Fixed parameters for all samples
CropLength=140
SampleSize=175000

# Input/output locations
RawDir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/rawData/raw_data_combined"
Outdir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/readProcessing/allDiploidFemales/out"
Adapters="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/readProcessing/adapters.fasta"

# Create outputdir if non-existent
mkdir -p "${Outdir}"

# iterate from samples list to creat running script
while IFS= read -r Base; do

    # Original input filesS
    Forward="${RawDir}/${Base}.R1.fq.gz"
    Reverse="${RawDir}/${Base}.R2.fq.gz"

    # Name of PBS script to be written
    Outfile="${Base}.readProcessing.pbs"

    (
    cat << EndOfJobScript
#!/bin/bash
#PBS -N ${Base}.psychids.readProcessingAllDiploidFemales.jun2026
#PBS -l select=1:ncpus=1:mem=36gb:scratch_local=600gb
#PBS -l walltime=39:59:00
#PBS -o /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/readProcessing/allDiploidFemales/out/${Base}.readProcessing.txt
#PBS -e /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/readProcessing/allDiploidFemales/err/${Base}.readProcessing.txt

set -euo pipefail

trap 'clean_scratch' TERM EXIT
cd "\${SCRATCHDIR}" || exit 1
export TMPDIR="\${SCRATCHDIR}"

# Copy input files to scratch using rsync
rsync -a "${Forward}" "\$SCRATCHDIR/"
rsync -a "${Reverse}" "\$SCRATCHDIR/"


# Create intermediate output folder
mkdir -p intermediateOutputs

# Define the local scratch copies explicitly
LocalForward="${Base}.R1.fq.gz"
LocalReverse="${Base}.R2.fq.gz"

# Define local seqtk binary
SEQTK="/storage/plzen1/home/pedroribeiro/Software/seqtk/seqtk"

##### Trim and filter #####

module add python27-modules/gcc
/software/python27-modules/software/python-2.7.6/gcc/bin/cutadapt \
-j 0 \
-b "file:${Adapters}" \
-B "file:${Adapters}" \
--nextseq-trim=10 \
--length=${CropLength} \
--minimum-length=${CropLength} \
-o intermediateOutputs/${Base}_filt_1.fastq.gz \
-p intermediateOutputs/${Base}_filt_2.fastq.gz \
"\${LocalForward}" "\${LocalReverse}"
module unload python27-modules/gcc

echo "cutadapt done" >&2

##### Sample reads #####

echo "using local seqtk" >&2

zcat "intermediateOutputs/${Base}_filt_1.fastq.gz" | \
sed -E "s/^@.*[ /]([12]):.*$/@${Base}\/\1/" | "\${SEQTK}" sample -s 10 - "${SampleSize}" | gzip -c > "intermediateOutputs/${Base}_sample_1.fastq.gz"

echo "zcat 1 done" >&2

zcat "intermediateOutputs/${Base}_filt_2.fastq.gz" | \
sed -E "s/^@.*[ /]([12]):.*$/@${Base}\/\1/" | "\${SEQTK}" sample -s 10 - "${SampleSize}" | gzip -c > "intermediateOutputs/${Base}_sample_2.fastq.gz"

echo "zcat 2 done" >&2

"\${SEQTK}" mergepe "intermediateOutputs/${Base}_sample_1.fastq.gz" "intermediateOutputs/${Base}_sample_2.fastq.gz" > "intermediateOutputs/${Base}_sample_12.fastq.gz"

echo "mergepe done" >&2

##### converts to FASTA #####

"\${SEQTK}" seq -A "intermediateOutputs/${Base}_sample_12.fastq.gz" > "${Base}_12.fasta"

echo "fasta conversion done" >&2

rm "${Base}.R1.fq.gz"
rm "${Base}.R2.fq.gz"

echo "raw reads were removed" >&2

# Copy outputs back
rsync -a "\$SCRATCHDIR/" "${Outdir}/" || export CLEAN_SCRATCH=false
EndOfJobScript
    ) > "${Outfile}"

    chmod +x "${Outfile}"
    echo "Script written: ${Outfile}"

done < samplesList.txt

SubmissionScript="submit_readProcessing_jobs.sh"

(
cat << 'EndOfSubmissionScript'
#!/bin/bash
EndOfSubmissionScript

while IFS= read -r Base; do
    echo "qsub ${Base}.readProcessing.pbs"
done < samplesList.txt

) > "${SubmissionScript}"

chmod +x "${SubmissionScript}"
echo "Submission script written: ${SubmissionScript}"
