#PBS -N psychids.allDiploidFemalesReMappingRefBuild.sep2026
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=50gb
#PBS -l walltime=00:30:00
#PBS -o /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/allDiploidFemales/out/allDiploidFemalesReMappingRefBuild.txt
#PBS -e /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/allDiploidFemales/err/allDiploidFemalesReMappingRefBuild.txt

set -euo pipefail

trap 'clean_scratch' TERM EXIT

##### system specific definitions #####
cd $SCRATCHDIR || exit 1
export TMPDIR=$SCRATCHDIR

Ref="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/REclustering/allDiploidFemales/out/PSYallDiploidFemalesConcatenated/output/seqclust/clustering/clusters/mappingRef/allDiploidFemalesConcatenatedRef.fasta"
OutRef="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/mapping/allDiploidFemales/out/refBuild"

mkdir -p "${OutRef}"

rsync -a "${Ref}" "$SCRATCHDIR/"

# Build bowtie2 index
module add bowtie2
bowtie2-build allDiploidFemalesConcatenatedRef.fasta allDiploidFemalesConcatenatedRef

rsync -a "$SCRATCHDIR/" "${OutRef}/"
