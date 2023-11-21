
# A07a5_samtools_stats.py ======================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_samstats(filepath):

    term_dict = {
        'raw total sequences': 'FilteredSeqCount',
        'bases mapped' : 'BasesMapped',
        'error rate': 'ErrorRate'
        }

    with open(filepath) as report:
        report_dict = {}
        for line in report:
            try:
                lhs, rhs = line.split(':')
            except ValueError:
                continue
            try:
                report_dict[term_dict[lhs]] = rhs.strip().split('\t')[0]
            except KeyError:
                pass
            
    return(report_dict)



# gather metadata --------------------------------------------------------------

metadata_well = pd.read_csv(filepath_wellmetadat)

filelist = metadata_well['A04c_txt_samstats']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_samstats = [parse_samstats(f) for f in filelist[boolean_fileexists]]
df_samstats = pd.DataFrame(list_samstats,
                        index=metadata_well['wellprefix'][boolean_fileexists])



# print percent files missing
print("number of target files: " + str(len(filelist)))
print("percent files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# column QC
print("Number of NAs per column:")
print(df_samstats.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_samstats.index.duplicated().sum()
print(ndupe)



# final export
df_samstats.to_csv("Metadata/A07a5_samstats.tsv", sep='\t')
print("done.\n\n")
