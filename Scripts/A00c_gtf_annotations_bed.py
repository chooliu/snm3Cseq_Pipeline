
# ==============================================================================
# A00c_gtf_annotations_bed.py 
# exports four annotation-related files to $ref_dir (reference genome)
# for mcds creation & down-stream analysis
# ==============================================================================

# recommend running interactively in python/Jupyter to check outputs
# works for .gtfs from GENCODE, but double check fmt for other sources

# note: if running interactively, you can source/export "snm3C_parameters.env"
# or, manually spec os.environ features by commenting out below 3 lines e.g.,
# os.environ['ref_dir'] = "/u/project/cluo/chliu/Genomes/IGVF_hg38_pluslambda"
# os.environ['ref_gtf'] = "/u/project/cluo/chliu/Genomes/IGVF_hg38_pluslambda/IGVF_GRCh38_plus_lambda.fa"
# os.environ['ref_chromsizes'] = "/u/project/cluo/chliu/Genomes/IGVF_hg38_pluslambda/chromsizes.tsv"



# load packages ----------------------------------------------------------------

import pandas as pd
import os



# load reference info ----------------------------------------------------------

os.chdir(os.environ['ref_dir'])
os.makedirs("annotations/", exist_ok=True)

gtf_file = pd.read_csv(os.environ['ref_gtf'],
                       comment = "#", delimiter="\t", header = None)
chrom_sizes = pd.read_csv(os.environ['ref_chromsizes'], sep = "\t", header = None)
chrom_sizes.columns = ['#chr', 'chrlen'] 



# genebody ---------------------------------------------------------------------
# .gtf to .bed (1-based --> 0 based)
# columns: chr, start, end, ENSG identifier
# note that this contains mitochondrial, ribosomal, lncRNAs, etc;
# this may or may not be desireable in downstream analyses

genebody = gtf_file[gtf_file.iloc[:, 2] == 'gene'].iloc[:, [0, 3, 4, 8]]
genebody.iloc[:, 1] = genebody.iloc[:, 1] - 1 # start changes to 0-pos
genebody.columns = ['#chr', 'start', 'end', 'annot']
genebody.reset_index(inplace=True, drop=True)

# extract info from the annot column (;, ")
genebody['gene'] = genebody['annot'].transform(lambda x: str(x).split('\"')[1])
genebody['symbol'] = genebody['annot'].transform(lambda x: str(x).split('\"')[5])
genebody['type'] = genebody['annot'].transform(lambda x: str(x).split('\"')[3])

# .gtf checks: should be zero
if sum(genebody.gene.duplicated()) != 0 | sum(genebody.start >= genebody.end) != 0:
  print("WARNING: check .bed outputs; was gene info was processed correctly from .gtf?")

# export ENSG --> Symbol
genebody.drop('annot', axis = 1, inplace = True)
genebody.to_csv("annotations/ensembl_to_symbol.tsv", sep = "\t", index = False)

# bed4 format for .allcools
# (drops strand info, which can be useful for other software)
genebody = genebody.iloc[:, 0:4]
genebody.to_csv("annotations/genebody.bed", sep = "\t", index = False)



# gene +/- 2kb -----------------------------------------------------------------
# above, but padding a 2kb region
# past manuscripts do this for higher mC modality coverage, but less interpretable

g2k = genebody.copy()
g2k.iloc[:, 1] = g2k.iloc[:, 1] - 2000
g2k.iloc[:, 2] = g2k.iloc[:, 2] + 2000

# but check +/-2kb still within chromosome length
# (low # genes affected, but may cause downstream issues)
g2k.loc[g2k.start < 0, 'start'] = 0

g2k = pd.merge(g2k, chrom_sizes, on = '#chr') 
filter_chrlen = g2k.end > g2k.chrlen
g2k.loc[filter_chrlen, 'end'] = g2k.chrlen[filter_chrlen]
g2k.drop('chrlen', axis = 1, inplace=True)
g2k.to_csv("annotations/geneslop2k.bed", sep = "\t", index = False)



# exonic -----------------------------------------------------------------------
# exon-level annotations, useful for STAR exon-only quant, potential ASE

exon = gtf_file[gtf_file.iloc[:, 2] == 'exon'].iloc[:, [0, 3, 4, 8]]
exon.iloc[:, 1] = exon.iloc[:, 1] - 1 # start changes to 0-pos
exon.columns = ['#chr', 'start', 'end', 'annot']

exon['gene'] = exon['annot'].transform(lambda x: str(x).split('\"')[1])
exon['transcript'] =  exon['annot'].transform(lambda x: str(x).split('\"')[3])
exon = exon.drop('annot', axis = 1).reset_index(drop=True)
exon.to_csv("annotations/exon.bed", sep = "\t", index = False)



# rRNA genes -------------------------------------------------------------------
# for some QC metrics after RNA alignments

rRNA = genebody.loc[gtf_file.iloc[:, 8].str.contains("rRNA"), :]
rRNA.to_csv("annotations/rRNA.bed", sep = "\t", index = False)


