# Repeatome comparisons in *Dahlica triquetrella* (RepeatExplorer2 + DESeq2)

Comparative repeat-content analysis of short-read libraries of the bagworm moth *Dahlica triquetrella* (Psychidae). Repeats are discovered de novo with RepeatExplorer2 (RE2), every sample is mapped back onto the RE2 clusters, and per-cluster abundance is compared between groups with DESeq2.

The repository is a record of the analysis as it was run, organised so that someone else can rerun it from raw reads or pick it up from any step.

## Design

```
raw reads ──► 1 trim + subsample ──► 2 concatenate ──► 3 RE2 clustering ──► 4 mapping reference
                    │                                                              │
                    └──────────────► 6 map every sample ◄── 5 bowtie2 index ◄──────┘
                                            │
                                            ▼
                              7 cluster x sample matrix ──► R: DESeq2
```

1. **Read processing.** Adapter trimming and cropping to exactly 140 bp (cutadapt), random subsampling of 175,000 reads per mate (seqtk, seed 10), interleaving and conversion to FASTA. Read names are rewritten to `@PSY_xxxx/1` and `@PSY_xxxx/2`.
2. **Concatenation.** All per-sample FASTAs of a comparison go into one file. Sample identity survives because of the read-name prefix from step 1.
3. **Clustering.** RE2 `seqclust` on the single concatenated FASTA. RE2 recovers each read's sample from the first 6 characters of its name (`--prefix_length 6`), so one clustering run produces one shared set of repeat clusters for all samples. This is why the clustering list (`*Concatenated.txt`) contains a single name.
4. **Mapping reference.** Contigs and TAREAN consensus sequences of every cluster are collected, and every header is prefixed with its cluster ID (`CL0001_…`). This prefix is what step 7 uses to collapse coverage per cluster.
5. **Index.** `bowtie2-build` on the mapping reference.
6. **Mapping.** Each sample is mapped with bowtie2, and per-base depth is computed with `bedtools genomecov`.
7. **Matrix.** Per-base depth is summed per cluster into a cluster x sample matrix, which is the DESeq2 input.

### Fixed parameters (identical in every comparison)

| Step | Setting |
|---|---|
| cutadapt | adapters on both mates (`-b`/`-B`), `--nextseq-trim=10`, `--length 140 --minimum-length 140` |
| subsampling | `seqtk sample -s 10`, 175,000 reads per mate |
| RE2 `seqclust` | `--paired --automatic_filtering --mincl 0.005 --prefix_length 6 --assembly_min 5 --taxon METAZOA3.0 --max_memory 1000000000 --cleanup` |
| bowtie2 | `-N 1 -X 101 --threads 2` |
| coverage | `bedtools genomecov`, per base |

## The three comparisons

| Folder | Question | Samples | Steps run | Reference |
|---|---|---|---|---|
| `comparison_males_full_re2_deseq2_pipeline` (key `geo`) | Do sexual males differ in repeat content between Swiss populations and Linz (Austria)? | 18 | 1–7 | its own |
| `comparison_mVSf_sexuals_full_re2_deseq2_pipeline` (key `mVSf_sexuals`) | Do sexual females follow the same geographic pattern as sexual males? | 29 | 1, 6, 7 | geo's reference |
| `comparison_diploid_females_full_re2_deseq2_pipeline` (key `allDiploidFemales`) | Does repeat content differ between sexual and parthenogenetic diploid females? | 18 (17 analysed) | 1–7 | its own |

The key in brackets is the name used in data paths and file names (`readProcessing/geo/`, `geo_cluster_matrix.tsv`, …).

**Why mVSf has no clustering.** It reuses the geo (male) reference, so that females are measured on the same cluster set as the males they are compared with. A reference does not need to be built from the samples mapped onto it. The two other comparisons have their own references by design: a reference built from one comparison's samples would bias another comparison's contrast.

## Repository layout

```
comparison_<name>/
├── hpc/   MetaCentrum (PBS) scripts, sample lists
└── R/
    ├── inputs/    cluster matrix (from step 7) + sample metadata
    ├── scripts/   one DESeq2 script
    └── outputs/   every table and figure the script produces
```

### Scripts per step

| Step | males (geo) | mVSf_sexuals | allDiploidFemales |
|---|---|---|---|
| 1 | `wrapperReadProcessing.sh` | `wrapperReadProcessing.sh` | `wrapperReadProcessing.sh` |
| 2 | `concatenateForClustering.sh` | – | `concatenateForClustering.sh` |
| 3 | `wrapperREclustering.sh` | – | `wrapperREclusteringDipFem.sh` |
| 4 | `psychids.geoCreateMappingRef.apr2026.sh` | – | `psychids.allDiploidFemalesCreateMappingRef.sep2026.sh` |
| 5 | `psychids.geoReMappingRefBuild.apr2026.pbs` | – | `psychids.allDiploidFemalesReMappingRefBuild.sep2026.sh` |
| 6 | `wrapperReMapping.sh` | `wrapperReMapping.sh` | `wrapperMappingAllDiploidFemales.sh` |
| 7 | `psychids.summarize_forDESeq.may2026.sh` | `psychids.summarize_forDESeq.may2026.sh` | `psychids.summarize_forDESeq.sep2026.sh` |

`PSYgeoConcatenated_clustering.pbs` (males) is the clustering job that `wrapperREclustering.sh` generated for the real run, kept as a record of the exact command and resources.

## Running the HPC part (MetaCentrum)

- **Wrappers** (`wrapper*.sh`) are run on a login node from inside `hpc/` (`bash wrapperX.sh`). They read `samplesList.txt` (or the one-line `*Concatenated.txt` for clustering) and write one PBS job per sample, plus a submission script.
- **Standalone scripts** (steps 4, 5 and 7) are submitted directly with `qsub`.
- **Step 2** is plain bash on the login node. It refuses to overwrite an existing concatenated FASTA.
- The directories named in `#PBS -o/-e` must exist before `qsub`, because PBS opens them before the script body runs.
- Run the steps in order. Each step reads the previous step's output from the data root.
- **Trimmed reads for mapping (males and mVSf).** Step 1 writes the trimmed reads (`*_filt_1/2.fastq.gz`) into `intermediateOutputs/`, but the males and mVSf mapping wrappers read them from `intermediateOutputs/filtered/`, where they were moved by hand. For a rerun, move them there or change `ReadsDir` in those wrappers.
- **Step 4 writes into the RE2 output.** It creates `contigsDir/`, `tareanDir/` and `mappingRef/` inside the `clusters/` folder it reads from, so rerunning it overwrites them.

### Paths

All scripts use absolute paths from the original account, which need to point to the new locations before a rerun:

- **Data root:** `/storage/brno12-cerit/home/pedroribeiro/PostPhd/projects/psychids` holds the reads and all outputs. Under it: `rawData/raw_data_combined/` (raw reads), and `repeatExplorer/{readProcessing,REclustering,mapping,prepareForDESeq2}/<key>/out/`.
- **Log root:** `/storage/plzen1/home/pedroribeiro/postPhd/projects/psychids` is where the PBS and live logs are written.

### Software

- **cutadapt:** loaded through the MetaCentrum module `python27-modules/gcc`, called by its full path in `wrapperReadProcessing.sh`.
- **seqtk:** must be installed by the user. Change the `SEQTK=` line in the three `wrapperReadProcessing.sh` files to point to your binary.
- **RepeatExplorer2:** Singularity image `repex_tarean_0.3.12-7a7dc9e.sif`. The clustering scripts copy it into scratch before running (see the copy line in each clustering script for the source path). MetaCentrum also provides RE2 images under `/cvmfs/singularity.metacentrum.cz/RepeatExplorer/`, which is the commented `SING_IMAGE` line in the scripts. Keep the version fixed if results are to be comparable.
- **bowtie2, samtools, bedtools:** taken from the MetaCentrum software stack. See the module lines in the step 5 and step 6 scripts for the exact versions.

### Resources

Clustering is the heavy step. The heaviest run was the diploid-female clustering; the males' clustering needed far less. The requested resources are in the `#PBS -l` lines of each clustering script. Every sample contributes the same number of reads (175,000 per mate), so the resources needed depend only on the number of samples.

## Running the R part

R ≥ 4.5 with DESeq2 (Bioconductor: `BiocManager::install("DESeq2")`). Each script installs any missing CRAN packages it needs.

1. Copy the comparison's matrix (`prepareForDESeq2/<key>/out/<key>_cluster_matrix.tsv`) into `R/inputs/`. The matrices used for the current results are already there.
2. In RStudio, set the working directory to that comparison's `R/` folder (Session › Set Working Directory › Choose Directory).
3. Open the script in `R/scripts/` and press **Source**. The script stops with a clear message if it is run from the wrong folder.

Every table and figure goes to `R/outputs/`, together with the fitted DESeq2 object and results (`*_dds_res.rds`) and `sessionInfo.txt` (R and package versions).

### Statistical choices (shared by all three scripts)

- **Counts.** The step 7 matrix holds summed per-base coverage, not read counts. Every read is exactly 140 bp and aligned end-to-end, so coverage ÷ 140 equals the number of mapped reads. The scripts convert to read counts before fitting, because DESeq2's negative-binomial model assumes counts of reads. Reads are counted per mate, so a pair whose two mates land in the same cluster counts twice.
- **Prevalence filter.** A cluster is tested only if it has at least 10 reads in at least N samples, where N is the size of the smallest group in that test (DESeq2 vignette recommendation). This removes clusters carried by one or two libraries, which cannot inform a group contrast but which DESeq2 can otherwise call significant.
- **Significance.** Benjamini–Hochberg adjusted p < 0.05. DESeq2's independent filtering runs at its default (alpha = 0.1).
- **Transform for PCA and QC.** `rlog`.
- **Library size.** Handled by DESeq2 size factors.

### Per comparison

**males (`geo`).** Design `~ geo_group`, Linz as reference (Swiss vs Linz). Two levels of analysis:

- All Swiss males against Linz.
- Each Swiss population (Nods, SwGuggis, RehetMid) against Linz separately.

From the per-population runs, the script lists the clusters significant in all three populations (`results_combined_geo.tsv`) and in at least two (`results_shared2_geo.tsv`), with a Venn diagram. The smallest group is Linz (n = 3), so N = 3 throughout.

**mVSf_sexuals.** The primary test is on females only: design `~ region`, Austrian as reference (Swiss vs Austrian). Sex is constant within the subset, so it is not in the model. The question is whether sexual females reproduce the geographic signal of sexual males, which validates using females as a baseline. A descriptive PCA of all 29 samples (males and females) is coloured by sex, population and region. That PCA has no model (`~ 1`) and uses only the rlog transform. The metadata lists PSY_0062, which has no column in the matrix. The script reports it and ignores it.

**allDiploidFemales.** Design `~ sex`, sexual as reference (parthenogenetic vs sexual; the data codes parthenogenetic as `asexual`). The comparison has 4 sexual females (PSY_0012, 0042, 0107, 0136) and 14 parthenogenetic females. PSY_0145 is excluded before fitting (`dropSamples`) because it belongs to another genus (*Pseudobankesia*), leaving 4 vs 13 and N = 4. This script also writes a dispersion plot and a sample-correlation heatmap as QC.

### Current results

| Comparison | Test | Clusters tested | Significant (padj < 0.05) |
|---|---|---|---|
| males | all Swiss vs Linz | 568 | 248 (44%) |
| males | Nods vs Linz | 568 | 258 |
| males | SwGuggis vs Linz | 568 | 238 |
| males | RehetMid vs Linz | 567 | 208 |
| mVSf_sexuals | females, Swiss vs Austrian | 567 | 252 (44%) |
| allDiploidFemales | parthenogenetic vs sexual | 734 | 3 (0.4%), all higher in parthenogenetic females |

Geography structures the repeatome strongly in both sexes. Reproductive mode barely does. The three significant clusters in diploid females are CL340, CL205 and CL280. Cluster IDs are only meaningful within one reference: CL340 in the diploid-female reference is unrelated to CL340 in the geo reference.
