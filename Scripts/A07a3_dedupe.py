
# A07a3_dedupe.py ==============================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd
import numpy as np

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

nulltable = np.array([pd.NA, pd.NA, pd.NA]) 

def parse_dedupe(filepath):
    try:
        data_dedupe = pd.read_csv(filepath, delimiter = "\t",
                         comment = "#", nrows = 1)[[
                             'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED', 'PERCENT_DUPLICATION'
                         ]].transpose()[0]
        return(data_dedupe)
    except:
        return(nulltable)

tidy_name_dict = {'PERCENT_DUPLICATION' : 'picard_perc_dupe',
                  'READ_PAIRS_EXAMINED' : 'picard_npairsin',
                  'UNPAIRED_READS_EXAMINED' : 'picard_nreadsin'}



# gather metadata --------------------------------------------------------------

filelist=metadata_well['A04a_log_picard']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picard = [parse_dedupe(f) for f in filelist[boolean_fileexists]]

df_picard = pd.DataFrame(list_picard,
                            index = metadata_well['wellprefix'][boolean_fileexists]
                           ).rename(columns = tidy_name_dict
                           ).drop("picard_npairsin", axis = 1)


# print percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# column QC
print("Number of NAs per column:")
print(df_picard.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_picard.index.duplicated().sum()
print(ndupe)



# final export
df_picard.to_csv("Metadata/A07a3_dedupe.tsv", sep='\t')
print("done.\n\n")
