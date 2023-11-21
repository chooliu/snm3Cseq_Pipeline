
# A07a7_contacts.py ============================================================
# setup ------------------------------------------------------------------------

import os
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_pairs(filepath, prefix = ""):
    return(pd.read_csv(filepath, delimiter="\t"))



# gather metadata --------------------------------------------------------------

filelist=metadata_well['A06a_3c_metadat']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_contacts = [parse_pairs(file) for file in filelist[boolean_fileexists]]
df_contacts = pd.concat(list_contacts).set_index(metadata_well['wellprefix'][boolean_fileexists])



# print percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# column QC
print("Number of NAs per column:")
print(df_contacts.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_contacts.index.duplicated().sum()
print(ndupe)



# final export
df_contacts.to_csv("Metadata/A07a7_contact_metadat.tsv", sep='\t')
print("done.\n\n")
