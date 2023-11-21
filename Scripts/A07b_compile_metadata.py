
# A07b_compile_metadata.py =====================================================
# assumes no changes to script output names from A07a*s

# setup ------------------------------------------------------------------------

import pandas as pd
    


# load tables ------------------------------------------------------------------
# basic aggregation based on wellprefix (should be index of each A07a* output)

def read_tbl_wrapper(filepath, prefix = ""):
    return(pd.read_csv(filepath, delimiter = "\t", index_col = 0).add_prefix(prefix))

target_files = ["Metadata/A07a1_trimming.tsv",
                "Metadata/A07a2_mappingrate.tsv",
                "Metadata/A07a3_dedupe.tsv",
                "Metadata/A07a4_global_mC_fracs.tsv",
                "Metadata/A07a5_samstats.tsv",
                "Metadata/A07a6_DNA_cov_chrXdivY.tsv", "Metadata/A07a6_DNA_cov_percent1x.tsv",
                "Metadata/A07a7_contact_metadat.tsv"]
                
# Note: joining all tables may yield
# "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
# which indicates that there's a metadata file containing duplicated wells
# this usually is a problem with script A07a4
metadata_mC = pd.concat([read_tbl_wrapper(f) for f in target_files], axis = 1)



# column QC
print("Number of NAs per column:")
print("(note, these may vary from A07a* scripts b/c those skip missing files)")
print(metadata_mC.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=metadata_mC.index.duplicated().sum()
print(ndupe)



# final export
metadata_mC.to_csv("Metadata/A07b_compiled_metadata.tsv", sep = "\t")
print("done.\n\n")
