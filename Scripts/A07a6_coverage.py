
# A07a6_coverage.py ============================================================
# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd
import numpy as np

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)



# extract autosomal ------------------------------------------------------------

target_chroms = ["chr" + str(i) for i in range(1, 99)]
autosomal_chroms = \
    pd.read_csv(os.environ['ref_chromsizes'],
                sep = "\t", header = None, index_col = 0)
autosomal_chroms = autosomal_chroms[autosomal_chroms.index.isin(target_chroms)]
total_autosomal_bases = autosomal_chroms.sum()
target_chroms=autosomal_chroms.index



# gather metadata: base-lvl unique coverage levels -----------------------------

print("processing autosomal num sites with at least 1-fold coverage.")
print("if any filenames printed below, potentially corrupt files:")
def parse_coverage_unique(filepath):
    try:
        percent_coverage = \
            pd.read_csv(filepath, delimiter = "\s+", header = None, index_col=1)
        percent_coverage = (percent_coverage.loc[target_chroms, 0].sum() / total_autosomal_bases).to_list()[0]
    except:
        print("potentially corrupt file:")
        print(filepath)
        percent_coverage = np.nan
    return(percent_coverage)

filelist=metadata_well['A04c_txt_covnsites']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_unique = [parse_coverage_unique(file) for file in filelist[boolean_fileexists]]
df_unique = pd.DataFrame(list_unique,
                        index = metadata_well['wellprefix'][boolean_fileexists])
df_unique.columns = ["CoveragePerc1x"]



# print percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# column QC
print("Number of NAs per column:")
print(df_unique.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_unique.index.duplicated().sum()
print(ndupe)



df_unique.to_csv("Metadata/A07a6_DNA_cov_percent1x.tsv", sep='\t')

print("done.\n\n")



# total coverage levels for chrX/chrY ------------------------------------------

print("processing total coverage levels per chrom.")
print("if any filenames printed below, potentially corrupt files:")
def parse_coverage_total(filepath):
    try:
        total_cov_by_chr = pd.read_csv(filepath, delimiter = "\s+", header = None, index_col=0)
        if any(total_cov_by_chr.index=="chrX") and (not any(total_cov_by_chr.index=="chrY")):
            coverage_XdivY = numpy.inf
        else:
            coverage_XdivY = total_cov_by_chr.loc['chrX', ] / total_cov_by_chr.loc['chrY', ]
            coverage_XdivY = coverage_XdivY.tolist()[0]
    except:
        print(filepath)
        coverage_XdivY = np.nan
    return(coverage_XdivY)

filelist=metadata_well['A04c_txt_covtot']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_total = [parse_coverage_total(file) for file in filelist[boolean_fileexists]]
df_total = pd.DataFrame(list_total,
             index = metadata_well['wellprefix'][boolean_fileexists])
df_total.columns = ["CoverageXdivY"]



# print percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# column QC
print("Number of NAs per column:")
print(df_total.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_total.index.duplicated().sum()
print(ndupe)



# final export
df_total.to_csv("Metadata/A07a6_DNA_cov_chrXdivY.tsv", sep='\t')
print("done.\n\n")
