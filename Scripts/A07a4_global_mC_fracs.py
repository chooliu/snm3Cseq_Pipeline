
# A07a4_global_mC_fracs.py =====================================================

# setup ------------------------------------------------------------------------

import glob
import pandas as pd
import os

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

attempt_dedupe = True # <-- attempt index repair?

# gather metadata --------------------------------------------------------------

filelist=[ "Metadata/A04d_mCfrac_" + str(batch_num) + ".tsv"
           for batch_num in pd.unique(metadata_well['batchnum']) ]
boolean_fileexists = [os.path.exists(f) for f in filelist]

list_mCfracs = [ pd.read_csv(f, delimiter = "\t") for f in pd.Series(filelist)[boolean_fileexists]] 
df_mCfracs = pd.concat(list_mCfracs)
df_mCfracs = df_mCfracs.rename(columns = {"Well" : "wellprefix"})



# print number files
print("number of target files: " + str(len(filelist)))

# column QC
print("Number of NAs per column:")
print(df_mCfracs.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_mCfracs.wellprefix.duplicated().sum()
print(ndupe)



# check for dupe wells ---------------------------------------------------------
# the A04d scripts can generate duplicated wells due to the "append" option
# or early terminated jobs which may result in NAs
# thus we can sort by number of CHs observed in desc order --> remove first dupes
if attempt_dedupe and ndupe != 0:
    print("attempting to dedupe...")
    df_mCfracs = df_mCfracs.sort_values('CH')
    df_mCfracs = df_mCfracs[~df_mCfracs.wellprefix.duplicated(keep = 'first')]



# final export
df_mCfracs = df_mCfracs.sort_values('wellprefix').reset_index(drop=True)
df_mCfracs.to_csv("Metadata/A07a4_global_mC_fracs.tsv", sep='\t', index = False)
print("done.\n\n")
