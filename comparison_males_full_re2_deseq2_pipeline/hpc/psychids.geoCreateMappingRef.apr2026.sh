#!/bin/bash
#PBS -N psychids.geoCreateMappingRef.apr2026
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=23:00:00
#PBS -o /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/geo/out/geoCreatMappingRef.txt
#PBS -e /storage/plzen1/home/pedroribeiro/postPhd/projects/psychids/repeatExplorer/mapping/geo/err/geoCreatMappingRef.txt

CLdirs="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/REclustering/geo/out/PSYgeoConcatenated/output/seqclust/clustering/clusters"
clusterList="/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids/repeatExplorer/REclustering/geo/out/PSYgeoConcatenated/output/seqclust/clustering/clusters/clusters.txt"

contigsDir="${CLdirs}/contigsDir"
tareanDir="${CLdirs}/tareanDir"
mappingRef="${CLdirs}/mappingRef"

mkdir -p "${contigsDir}"
mkdir -p "${tareanDir}"
mkdir -p "${mappingRef}"

##### copy contigs.fasta and consensus_dimer.fasta #####
while IFS= read -r clusterID; do

    clusterDir="dir_${clusterID}"

    if [ -f "${CLdirs}/${clusterDir}/contigs.fasta" ]; then
        cp "${CLdirs}/${clusterDir}/contigs.fasta" "${contigsDir}/${clusterID}.contigs.fasta"
        echo "Copied contigs for ${clusterID}" >&2
    else
        echo "Missing contigs.fasta for ${clusterID}" >&2
    fi

    if [ -f "${CLdirs}/${clusterDir}/tarean/consensus_dimer.fasta" ]; then
        cp "${CLdirs}/${clusterDir}/tarean/consensus_dimer.fasta" "${tareanDir}/${clusterID}.consensus_dimer.fasta"
        echo "Copied TAREAN consensus for ${clusterID}" >&2
    fi

done < "${clusterList}"

##### edit contigs headers #####
for fasta in "${contigsDir}"/*.fasta; do
    [ -e "${fasta}" ] || continue
    sed -i 's/^>\(CL[0-9]\+\)\(Contig[0-9]\+\)$/>\1_\2/' "${fasta}"
    echo "Edited contig headers in $(basename "${fasta}")" >&2
done

##### edit TAREAN headers #####
for fasta in "${tareanDir}"/*.fasta; do
    [ -e "${fasta}" ] || continue

    fileBase=$(basename "${fasta}" .consensus_dimer.fasta)
    clusterNum=$(echo "${fileBase}" | sed 's/^CL0*//')
    clusterShort="CL${clusterNum}"

    sed -i "s/^>/>${clusterShort}_/" "${fasta}"
    echo "Edited TAREAN headers in $(basename "${fasta}")" >&2
done

##### concatenate edited contigs #####
allContigs="${contigsDir}/geo.allContigs.fasta"

> "${allContigs}"
for fasta in "${contigsDir}"/CL*.contigs.fasta; do
    [ -e "${fasta}" ] || continue
    cat "${fasta}" >> "${allContigs}"
done

echo "Created combined contigs file: $(basename "${allContigs}")" >&2

##### concatenate edited TAREAN consensuses #####
allTarean="${tareanDir}/geo.allTarean.fasta"

> "${allTarean}"
for fasta in "${tareanDir}"/CL*.consensus_dimer.fasta; do
    [ -e "${fasta}" ] || continue
    cat "${fasta}" >> "${allTarean}"
done

echo "Created combined TAREAN file: $(basename "${allTarean}")" >&2

##### copy combined files to mappingRef #####
cp "${allContigs}" "${mappingRef}/"
cp "${allTarean}" "${mappingRef}/"

echo "Copied combined contigs and TAREAN files to mappingRef" >&2

##### create final mapping reference #####
cat "${mappingRef}/geo.allContigs.fasta" \
    "${mappingRef}/geo.allTarean.fasta" \
    > "${mappingRef}/geoConcatenatedRef.fasta"

echo "Created final mapping reference: geoConcatenatedRef.fasta" >&2
