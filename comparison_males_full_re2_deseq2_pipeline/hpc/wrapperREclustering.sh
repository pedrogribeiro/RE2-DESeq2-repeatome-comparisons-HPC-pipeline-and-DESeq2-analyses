#!/bin/bash

PrefixLength=6

# Loop through all fastas
while IFS= read -r Base; do

    # Define file paths
    Fasta="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/readProcessing/geo/out/${Base}.fasta"

    Outdir="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/REclustering/geo/out/${Base}"
    mkdir -p "${Outdir}"
    Outfile="${Base}_clustering.pbs"

    # Write the PBS script
    (
    cat << EndOfJobScript
#!/bin/bash
#PBS -N ${Base}.psychids.geoREclustering.apr2026
#PBS -l select=1:ncpus=8:mem=72gb:scratch_local=1400gb:cl_bee=True:singularity=True
#PBS -l walltime=23:59:00
#PBS -o /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/REclustering/geo/out/${Base}.REclustering.txt
#PBS -e /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/REclustering/geo/err/${Base}.REclustering.txt

trap 'clean_scratch' TERM EXIT

export TMPDIR=\${SCRATCHDIR}

rsync -a "${Fasta}" "\$SCRATCHDIR/"
cd "\${SCRATCHDIR}"

# define local fasta explicitly after copy to scratch
LocalFasta="${Base}.fasta"

# define variables
#SING_IMAGE=/cvmfs/singularity.metacentrum.cz/RepeatExplorer/repex_tarean.sif

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "\$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

#pull the latest RE image; check it at https://cloud.sylabs.io/library/repeatexplorer/default/repex_tarean
singularity pull --arch amd64 library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e

#set SINGULARITY variables for runtime data
export SINGULARITY_CACHEDIR=\${SCRATCHDIR}
export SINGULARITY_LOCALCACHEDIR=\${SCRATCHDIR}
export SINGULARITY_TMPDIR=\${SCRATCHDIR}

singularity exec --bind \${SCRATCHDIR}:/data/ --bind \$SCRATCHDIR:/tmp/ \${SCRATCHDIR}/repex_tarean_0.3.12-7a7dc9e.sif seqclust --paired --automatic_filtering --mincl 0.005 --cpu 8 --prefix_length ${PrefixLength} \
--output_dir=/data/output/ --max_memory 1000000000 --cleanup --assembly_min 5 --taxon METAZOA3.0 /data/\${LocalFasta}

rm -r \${SCRATCHDIR}/repex_tarean_0.3.12-7a7dc9e.sif
rsync -a \${SCRATCHDIR}/output "${Outdir}/" || export CLEAN_SCRATCH=false

EndOfJobScript
) > "${Outfile}"

    chmod +x "${Outfile}"
    echo "Script written: ${Outfile}"

done < geoConcatenated.txt
