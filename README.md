# sn-m3C-seq Bioinformatics Pipeline



![Overview of sn-m3C-seq pipeline.](./Documentation/snm3C_overview.png)

**sn-m3C-seq** simultaneously profiles DNA methylome (mC) and chromatin contacts within a single nucleus. As DNA methylation is unaltered during chromatin conformation capture (3C) protocols, we perform 3C/HiC-style _in situ_ reactions followed by nucleus isolation into 384-well plates, bisulfite conversion, and sequencing library generation on the 3C ligation products.  

**What's in this repo?**

* **An updated version of TAURUS-MH** to go from plate-level .fastqs &rarr; cell-level alignments (.bam), quality control metrics, and core mC & 3C features.
* For the **methylome**, the primary features are the methylation counts and coverage for each cytosine in the genome (**.allc** file), which can then be aggregated into methylation levels of different genomic intervals (e.g., 100kb-bins, genes; **.mcds** file).
* For **chromatin conformation**, the features are read-level contacts (similar to **.pairs**), which can be aggregated into counts in different pairs of bins (10kb, 25kb, 100kb aggregated counts in sparse contact matrices _a la_ **.cool** files).

**What _isn't_ in this repo?**

* Cell-level QC steps, which can be celltype- and study-specific. (Some suggestions on the [Detailed Overview](./Documentation/detailed_overview.md) page and within [past manuscripts](#technology-references).)
* Downstream analysis of mC and 3C following basic quantification (e.g., feature selection/calling, clustering, hypothesis testing, imputation).

**Quick Start**
* Clone (`git clone`) or download this repo as a .zip file [from the releases page](https://github.com/chooliu/snm3Cseq_Pipeline/releases).
* install dependencies in .yaml: `conda env create -f Documentation/snm3Cseq_taurus.yml`
* customize the genome alignment and input sample metadata (either using Jupyter notebooks `A00*` and `A01*` *OR* by editing Scripts `A00*` and `A01*`)
* run submission scripts (`.sub`) in order, paying attention to parallelization task IDs `-t 1-Nplates` and `-t 1-Nbatches` (24 wells/batch by default). For convenience, this can be done via the [Notebooks](./Notebooks), or the full list of submission commands is listed at [Documentation/submission_helper.txt](./Documentation/submission_helper.txt).
* read the [Detailed Overview](./Documentation/detailed_overview.md) for FAQs and common pitfalls.

---

## Links

### Related Pipelines 

* Existing quantification pipelines, mainly from our collaborators at/formerly at the Salk Institute:
  - [TAURUS-MH](https://github.com/dixonlab/Taurus-MH): read-splitting based mapping pipeline historically used on sn-m3C-seq &rarr; developed specifically for this assay primarily by Dr. Dongsung Lee.
  - [YAP (Yet Another Pipeline)](https://hq-1.gitbook.io/mc/): supports sn-m3C-seq and additional related assays (e.g., mC, mCT, mCAT-seq), including new protocol based on restriction-site splitting for m3C. Snakemake. Developed primarily by the Ecker Lab/Drs. Hanqing Liu & Jingtian Zhou.

* Downstream analysis tools:
  - [allcools](https://lhqing.github.io/ALLCools): Typically used by our group for mC downstream analysis. Helpful to review for .allc and .mcds descriptions.
  - [scHiCluster](https://github.com/zhoujt1994/scHiCluster): Contact pairs to binned contact count matrices at various resolution, including imputation. Our groups have also worked with the Ma Lab's [Higashi](https://github.com/ma-compbio/Higashi/) & [FastHigashi](https://github.com/ma-compbio/Fast-Higashi).

### Technology References

* Flagship assay paper and closely related methylome-only reaction:
  - [sn-m3C-seq](https://pubmed.ncbi.nlm.nih.gov/31501549/): Lee DS, et al. Simultaneous profiling of 3D genome structure and DNA methylation in single human cells. Nature methods. 2019 Oct;16(10):999-1006.
  - [snmC-seq2](https://pubmed.ncbi.nlm.nih.gov/30237449/): Luo, C. et al. Robust single-cell DNA methylome profiling with snmC-seq2. Nat. Commun. 9, 7–12 (2018). [Note: we are now on [snmC-seq3](https://www.protocols.io/view/snm3c-seq3-cwxuxfnw.html) but only a wet lab protocol citation exists.]

* Our library structure is described in the [Detailed Overview](./Documentation/detailed_overview.md) and a [seqspec](https://igvf.github.io/seqspec/specs/sn-m3C-seq/spec.html). 

* Additional examples of applications/publications using sn-m3C-seq:
  - [Liu, et al. (2023) preprint](https://www.biorxiv.org/content/10.1101/2023.04.16.536509v1): Liu H, et al. Single-cell DNA Methylome and 3D Multi-omic Atlas of the Adult Mouse Brain. bioRxiv. 2023 Apr 18.
  - [Heffel, et al. (2022) preprint](https://www.biorxiv.org/content/10.1101/2022.10.07.511350v1): Heffel MG, et al. Epigenomic and chromosomal architectural reconfiguration in developing human frontal cortex and hippocampus. bioRxiv. 2022:2022-10.

### Example Datasets

* The example in scripts / detailed overview is for a fibroblast &rarr; induced pluripotent stem cell (iPSC) differentiation experiment; the underlying data to soon be published for the IGVF Consortium (NHGRI).
* Some raw and processed datasets are linked/publicly available via our lab website [luogenomics.github.io](https://luogenomics.github.io/data/). For example, note the `_allc.tar.gz`  (cytosine-level) and `_contacts.tar.gz` (contact pairs) files for [GEO GSM6596812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6596812). 

### Acknowledgements

* Thanks to original TAURUS-MH workflow developers, including Drs. Dongsung Lee, Chongyuan Luo, and Jesse Dixon.
* Luo Lab collaborators/members for pipeline testing and feedback, namely Zach von Behren and Min Jen Tsai.
