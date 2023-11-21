
# A07a1_trimming.py ============================================================

# setup ------------------------------------------------------------------------

import re
import pandas as pd
import glob
import os

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_fastp_report(filepath):
    jsonfile = pd.read_json(filepath)
    dict_out = {
        'nreads_pretrim' : jsonfile['summary']['before_filtering']['total_reads'],
        'percreads_passtrim' : jsonfile['summary']['after_filtering']['total_reads'] /
              jsonfile['summary']['before_filtering']['total_reads'],
        'q20_pretrim' : jsonfile['summary']['before_filtering']['q30_rate'],
        'q20_posttrim' : jsonfile['summary']['after_filtering']['q30_rate'],
        'r1_len' : jsonfile['summary']['after_filtering']['read1_mean_length'],
        'r2_len' : jsonfile['summary']['after_filtering']['read2_mean_length'],
        'gc_perc' : jsonfile['summary']['after_filtering']['gc_content']}
    return(dict_out)



# gather metadata --------------------------------------------------------------

filelist=metadata_well['A03a_json_fastp']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fastp = [parse_fastp_report(f) for f in filelist[boolean_fileexists]]
df_fastp = pd.DataFrame(list_fastp,
                        index=metadata_well['wellprefix'][boolean_fileexists])

# print percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# column QC
print("Number of NAs per column:")
print(df_fastp.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_fastp.index.duplicated().sum()
print(ndupe)



# final export
df_fastp.to_csv("Metadata/A07a1_trimming.tsv", sep='\t')
print("done.\n\n")
