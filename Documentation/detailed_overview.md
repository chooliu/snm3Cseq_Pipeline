
# Detailed Overview

![Overview of sn-m3C-seq pipeline.](../Documentation/snm3C_overview.png)



## Table of Contents

- [A00_environment_and_genome_setup](#a00_environment_and_genome_setup)
- [A01_mergefastq_preptargets](#a01_mergefastq_preptargets)
- [A02_demultiplex](#a02_demultiplex)
- [A03_trim](#a03_trim)
- [A04_map_bismark](#a04_map_bismark)
- [A05_quantify_mC](#a05_quantify_mc)
- [A06_quantify_contacts](#a06_quantify_contacts)
- [A07_compile_metadata](#a07_compile_metadata)
- [Final Directory Structure](#final-directory-structure)
---

<br><br><br>
# A00_environment_and_genome_setup

* Clone (`git clone`) or download this repo as a .zip file [from the releases page](https://github.com/chooliu/snm3Cseq_Pipeline/releases). Rename the resulting folder to an informative "project directory" to hold scripts & outputs.

* Install the dependencies in `Documentation/snm3Cseq_taurus.yaml`. Conda installation highly recommended for convenience/reproducibility.
```
module load anaconda3 # or otherwise access conda
conda env create -f Documentation/snm3Cseq_taurus.yml
```

* Prep the following directory structure in a new "project directory" where output files will be stored.
```
# dir_proj=/u/project/cluo_scratch/chliu/IGVF_iPSC_snm3Cseq_YZCL47
# mkdir $dir_proj; cd $dir_proj
# mkdir Metadata Notebooks Scripts
```
* ⚠🚨 Using the included Juypter `Notebooks` or a text editor, customize the project-specific parameters in `snm3C_parameters.env`, submission scripts (`Scripts/*.qsub`), and `A01b.py` + `A01c.py`, especially checking for:
    - (i) compatibility with your compute/scheduler infrastructure and sequencing depth\
    (resource suggestions assume ⪅ 3 million total reads/well average; generally 32 plates or fewer/NovaSeq S4)
    - (ii) genome/reference organism
    - (iii) **job array ranges**, which vary based on the number of 384-well plates profiled\
    \
    The array range `-t 1-Nplates` is used for tasks performed at the plate level (e.g., demultiplexing), but `-t 1-Nbatches` for more intensive tasks like alignment. By default, each "batch" contains 24 cells (one wellplate row).

* Run customized `*.qsub` submission scripts `A00a`, `A00b`, ..., `A01a`, etc., starting with reference genome prep (`.fasta` sequence and `.gtf` annotations) with the following commands.


### A00 Commands to Run ⭐

Each section has a code block showing typical `qsub` commands. These can be run all at once if `-hold_jid`/`-hold_jid_ad` or equivalent is supported by the scheduler. 

```
# check desired ref genome, paths, & snm3C_parameters.env before running
qsub Scripts/A00a_genome_dl_index.sub  
qsub Scripts/A00b_genome_prep_bismark.sub
qsub Scripts/A00c_annotations_bed.sub
```

**Always submit from within the project directory.** (Or if `qsub -cwd` flag not supported, edit to use absolute filepath for `.env` file, followed by `cd $dir_proj`)

<br><br><br>
# A01_mergefastq_preptargets

### A01 Commands to Run ⭐

```
qsub Scripts/A01a_merge_lanes.sub # * -t 1-Nplates
qsub Scripts/A01b_plate_metadata.sub
```

* Our raw data is demultiplexed at the plate-level (i5/i7 index reads) and base-called via Basespace/`bcl2fastq`. The resulting nomenclature is effectively `Plate_L00*_R1.fq.gz` (Read 1) and `Plate_L00*_R2.fq.gz` (Read 2). The "Plate" name usually contains study-specific metadata separated by dashes that are useful to retain.
* **A01a:** We generally assume no substantial lane-effect and concatenate a given plate's files into `fastq_raw/Plate_R1.fq.gz` and `fastq_raw/Plate_R2.fq.gz`. (Make sure the Plate identifier is unique.)
* **A01b:** Extract plate-level metadata from `.fq` names.
* **A01c:** The anticipated filepath intermediate and final file outputs are predictable for each well, so we expand the plate-level metadata into well-level metadata. Each well corresponds to one nucleus/cell, and is uniquely identified by its "`wellprefix`" (Plate_Well, where Well is the 384-well plate position in {A1, A2, ..., A23, A24, B1, ..., P23, P24}).
* ⚠🚨️ **The number of wells per batch is also defined in A01c.** There may be a tradeoff between the number of tasks and resources per task needed depending on the architecture of your cluster (e.g., if getting `h_rt=24:00:00` in bismark mapping is bottleneck, can switch to e.g., 12 wells/batch & 12 hours with twice as many jobs in the job array.) Time and resource suggestions are given in each submission script/Notebook.

### Example of Plate & Well Files
```
# one plate:
20231005-3C29D12-Pos1-B06_S16_L001_R1_001.fastq.gz
20231005-3C29D12-Pos1-B06_S16_L001_R2_001.fastq.gz
...
20231005-3C29D12-Pos1-B06_S16_L008_R1_001.fastq.gz
20231005-3C29D12-Pos1-B06_S16_L008_R2_001.fastq.gz
```

In A01, these are concatenated across sequencing lanes into `20231005-3C29D12-Pos1-B06_S16_R1.fastq.gz` and `_R2.fastq.gz` (verify that this filename still distinguishes the plate from others).

Plate names are parsed into `Metadata/A01b_plate_metadata.csv` (the dash-separated line, timepoint info specific for our 32-plate IGVF pilot experiment):
```
        plate                            dateseq        sample    sort    plateindex    line    time    platenum
1    20231005-3C39D12-Pos2-A09_S8    20231005    3C39D12    Pos2    A09_S8    C39    12    1
2    20231005-3C39D5-Pos1-C05_S26    20231005    3C39D5    Pos1    C05_S26    C39    5    2
...
31    20231005-3C39D9-Pos1-C09_S30    20231005    3C39D9    Pos1    C09_S30    C39    9    31
32    20231005-3C39D9-Pos2-A07_S6    20231005    3C39D9    Pos2    A07_S6    C39    9    32
```
`Metadata/A01c_well_filepaths.csv` has one row for well, where `Nplates*384` total wells are expected. Each well is uniquely identified by it's "`wellprefix`". The first wellplate row of `20220831-IGVF-A10-D0-C02_S18` constitute the first job task `batchnum` by default (`Nbatches = Ncells/24`).

```
    wellprefix                        batchnum        04a_bam_final
1    20231005-3C29D1-Pos1-C04_S25_A1    1    mapping_bismark/20231005-3C29D1-Pos1-C04_S25_A...
2    20231005-3C29D1-Pos1-C04_S25_A2    1    mapping_bismark/20231005-3C29D1-Pos1-C04_S25_A...
...
12286    20231005-3C39D9-Pos2-A07_S6_P22    512    mapping_bismark/20231005-3C39D9-Pos2-A07_S6_P2...
12287    20231005-3C39D9-Pos2-A07_S6_P23    512    mapping_bismark/20231005-3C39D9-Pos2-A07_S6_P2...
12288    20231005-3C39D9-Pos2-A07_S6_P24    512    mapping_bismark/20231005-3C39D9-Pos2-A07_S6_P2...
```

### A01 Troubleshooting Notes
* Explicit filepath specification is used to make one of the most manually time-consuming challenges of not using a workflow manager--checking file outputs/partially failed jobs--slightly less painful.
    - If a batch job didn't finish for whatever reason (e.g., some wells have unusually high read depth; node crashes), code like the following can help identify the job array to re-run:
    `metadat_well <- read.tsv("Metadata/A01c_well_filepath.tsv"); metadat_well$batchnum[sapply(targets$A04a_bam_bismark_PE, file.size) > 0]` in R
    - or analogous commands like`[os.path.getsize(f) > 0 for f in list_of_filepaths]` in python, or  `[[ ! -s $filepath ]]` looping through filepath in shell (see script `A02b` for examples.)
* Where possible, steps are ideally performed in POSIX-independent command line, or if needed python (to maintain conda environment/version control). There's thus some clunky components of each script: namely, parsing the filepaths with `$basharrays{[@]}` and `query_metadat` function in most scripts (finds column name in `metadat_well` &rarr; returns target rows associated with the batchnum or platenum of interest).
* If a spike-in bisulfite conversion control is used, add to reference genome (e.g., Lambda phage sequence, see `A00a` script).
* A study could also have on-going data collection. Specifying a new A01b and A01c targets file in the `snm3Cseq_parameters.env` file helps organize processing of newly assayed plates; however, it just may be tidier to just start a new project directory &rarr; merge afterwards.
* Alternatives to manually editing each `.sub` script to modify `#$ -t 1-Nplates` or `1-Nbatches`:
    - One could do a sed replace like the following: `sed -i.orig '/^#/ s/1-512/1-256/' Scripts/*.sub` (changed `A01c` to run 512 wells/batch &rarr; 256 wells/batch and now changing job arrays to reflect this; modifying every submission script in place but saving original with `.orig`).
    - If the number of plates differs or just one batch failed, setting `-t` via command line will supersede the range in the script. e.g.,  `qsub -t 2 Scripts/A01a_merge_lanes.sub` (re-runs just plate 2 even though `-t 1-Nplate` is in the script).

---

<br><br><br>
# A02_demultiplex

![Layout of snm3C library sequences.](../Documentation/library_structure.png)
For detailed sequences, please see our library's [seqspec](https://igvf.github.io/seqspec/specs/sn-m3C-seq/spec.html).

### A02 Commands to Run ⭐

```
qsub Scripts/A02a_demultiplex_fastq.sub # * -t 1-Nplates
qsub Scripts/A02b_check_demultip.sub #
qsub Scripts/A02c_fastqc_demultip_fastq.sub
```

* Each Read 1 begins with a 8bp cell barcode specific to one of the 384 wells in the plate, and by extension one nucleus/cell (barring relatively rare doublets or empty wells from flow cytometry sorting; likely <3%).
* The sequences and their corresponding positions are predefined (`Scripts/A02a_cellbarcodes_subset[1,2].fa`) and repeated across plates.
* A02a: We check the 8bp barcode for each read-pair, demultiplex the plate-level `fastq_raw/Plate_R1.fastq.gz` &rarr; well-level `fastq_demultip/Plate_Well_indexed_R1.fastq.gz` and `fastq_demultip/Plate_Well_indexed_R2.fastq.gz` files, for the 384 wells in {A1, A2, ..., P23, P24}.
* A02b: Check if expected outputs present, calculate number of empty wells.
* A02c: Apply `FastQC` on random subset of wells. Due to bisulfite conversion/cell barcodes, take % CG metrics & base composition plots with a grain of salt.

### A02 Troubleshooting Notes
* Two text summaries `fastq_demultip/Plate_summary_[1,2].txt` are printed to show the number of cell barcodes detected for each plate. The number of "unassigned" barcodes should be around ~50% because the cell barcodes safelist is queryed in two subsets to decrease resource demands. We typically see 1-3% reads per plate total unassigned.
* Wells with no cell barcodes detected will have no .fastq.gz output; some (typically) small number of empty wells is OK.
* I call the 8bp a "cell barcode" to unify our assay with convention, but past snmC pipelines/papers may call it an "index" (hence the "`indexed`" in the demultiplexed filenames). 
* We currently only check of exact barcode matches; may accomodate a 1-2bp mismatch in the future.

---

<br><br><br>
# A03_trim


### A03 Commands to Run ⭐
```
qsub Scripts/A03a_trimming_fastp.sub # † -t 1-Nbatchnum
qsub Scripts/A03b_check_trimmed.sub
qsub Scripts/A03c_fastqc_trimmed.sub
```

* Now that the sequences are demultiplexed, we remove the cell barcode, random hexamer priming sequence, artifactual adaptase tail in our expected library structure (figure above). In addition, we want to trim adapter/the random 9H primer sequences (due to mean insert sizes ~250-300, potentially present in both 5'/3' ends), low complexity sequences (polyN), and Q-score < 20 regions.
* There are scenarios where Read 1 passes QC to some minimum length, but Read 2 does not. These "trimming singletons" have to be mapped in single-end mode in subsequent steps. Four files will thus be generated in `fastq_trimmed`: R1 with a mate, R2 with a mate, R1 trimming singletons, R2 trimming singletons.
* **A03a.** Trimming via `fastp`.
* **A03b.** Checks all expected output files present.
* **A03c.** Re-runs `FastQC` for comparison to pre-trimming sequences. Should see the majority of cells pass QC, and shifts towards higher mean Q-Score, lower % adapter, and more even base compositions along position in the read.


### A03 Troubleshooting Notes:
* The "H" sequence refers to A, C, and T nucleotides.
* The contaminating sequence list is in `Scripts/A03a_adapter_sequences.fa`. Illumina adapter contamination is most common (Smart-primers rare). I included reverse compliments in case of very small insert size libraries.
* `FastQC` occasionally reports wells with high Illumina adapter sequence (often >50% adapter &rarr; still ~20-30% post-trimming) that are not removed by other trimming programs (incl. cutadapt, TrimGalore) or multiple rounds of trimming (comparable residual % adapter, but much higher resource requirements). Ancedotally, these wells often seem to have some extremely small insert sizes (low/fragmented template and majority artifacts/adapter-dimer?), as well as fewer absolute numbers of demultiplexed reads (by 10 to 100-fold). I don't currently calculate adapter % as a QC metric as these wells usually thus get excluded in downstream QC steps by virtue of read count or mapping rates.


---

**(Writing-in progress for Steps A04-A07!)**
