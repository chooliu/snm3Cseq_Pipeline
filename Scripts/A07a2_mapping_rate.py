
# A07a2_mapping_rate.py ========================================================

# setup ------------------------------------------------------------------------

import os
import glob
import re
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)
    
def parse_bismark_report(filepath):

    """
    parse bismark.txt output
    adapted from YAP @ https://github.com/lhqing/cemba_data to include PE & SE output
    commented out term_dict lines of limited interest
    note that paired-end metrics usually yield fragments, versus reads
    """

    term_dict = {
        'Sequence pairs analysed in total': 'TotalReadPairsIn',
        'Sequences analysed in total': 'TotalReadsIn',
        'Number of paired-end alignments with a unique best hit': 'UniqueMappedPairs',
        'Number of alignments with a unique best hit from the different alignments': 'UniqueMappedReads',
        'Mapping efficiency': 'MappingRate',
        # # other potential metrics, not usually used
#         'Sequence pairs with no alignments under any condition': 'UnmappedPairs',
#         'Sequences with no alignments under any condition': 'UnmappedReads',
#         'Sequences did not map uniquely': 'AmbigReads',
#         'Sequence pairs did not map uniquely': 'AmbigPairs',
#         'CT/GA/CT': 'ReadsOT',
#         'GA/CT/CT': 'ReadsOB',
#         'GA/CT/GA': 'ReadsCTOT',
#         'CT/GA/GA': 'ReadsCTOB',
#         'CT/CT': 'ReadsOT',
#         'CT/GA': 'ReadsOB',
#         'GA/CT': 'ReadsCTOT',
#         'GA/GA': 'ReadsCTOB',
#         'Total number of C\'s analysed': 'TotalC',
#         'C methylated in CpG context': 'BismarkmCGRate',
#         'C methylated in CHG context': 'BismarkmCHGRate',
#         'C methylated in CHH context': 'BismarkmCHHRate',
#         'C methylated in unknown context (CN or CHN)' : 'BismarkmCNCHNRate',
#         'C methylated in Unknown context (CN or CHN)' : 'BismarkmCNCHNRate'
        }

    with open(filepath) as report:
        report_dict = {}
        for line in report:
            try:
                lhs, rhs = line.split(':')
            except ValueError:
                continue
            try:
                report_dict[term_dict[lhs]] = rhs.strip().split('\t')[0].strip('%')
            except KeyError:
                pass
            
    return(report_dict)



# prep file lists --------------------------------------------------------------

# extract bismark logs that exist
target_filepath_columns = metadata_well.columns[
    metadata_well.columns.str.contains("A04a_bismarktxt_")]
filelist = metadata_well.loc[:, target_filepath_columns].values.tolist()
filelist = [f for sublist in filelist for f in sublist]
boolean_fileexists = [os.path.exists(f) for f in filelist]

# print percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))

# well list
# anticipate 10 unique logs per well (2 full length, up to 6 splits)
welllist = []
for w in metadata_well.wellprefix:
    welllist.extend([w]*10)

# loop through files
print("reading bismark logs...")
list_bismarkmap = [parse_bismark_report(f) for f in pd.Series(filelist)[boolean_fileexists]]

df_mapping_byalign = pd.DataFrame(list_bismarkmap).apply(pd.to_numeric)
df_mapping_byalign.index=pd.Series(welllist)[boolean_fileexists]

# pairing metadata by alignment type

print("pairing metadata with split/alignment type...")
df_mapping_byalign['Alignment'] = pd.Series(
    target_filepath_columns.str.replace("A04a_bismarktxt_", "").to_list() \
    * metadata_well.shape[0])[boolean_fileexists].to_list()

alignnames_presplit = ["R1p", "R2p", "R1trims", "R2trims"]

# can optionally save this detailed breakdown by alignment source
# df_mapping_byalign.to_csv("Metadata/A07a2_mappingrate_detailed.tsv", sep='\t')



# final metadata out -----------------------------------------------------------
# join the TAURUS pre-split (full R1 and R2 mapping), then post-split

print("now grouping by wellprefix...")

df_final = pd.DataFrame(index = metadata_well.wellprefix)

df_presplit = \
    df_mapping_byalign[df_mapping_byalign['Alignment'].isin(alignnames_presplit)
                      ].reset_index().groupby('index').agg('sum')
df_final['NumReadsIn'] = df_presplit.TotalReadsIn
df_final['UniqueMappedReads_PreSplit'] = df_presplit.UniqueMappedReads
df_final['MappingRate_PreSplit'] = df_presplit.UniqueMappedReads/df_presplit.TotalReadsIn

df_postsplit = \
    df_mapping_byalign[~df_mapping_byalign['Alignment'].isin(alignnames_presplit)
                      ].reset_index().groupby('index').agg('sum')

df_final['NumReadsIn_PostSplit'] = df_postsplit.TotalReadsIn
df_final['UniqueMappedReads_PostSplit'] = df_postsplit.UniqueMappedReads
df_final['MappingRate_PostSplit'] =  df_postsplit.UniqueMappedReads/df_postsplit.TotalReadsIn

# final combined mapping rates
# Nmappre + Nmappost/Nsplits; Nsplits = (Ninpost/Nunmapped) ~ 2-3 per read

print("doing final summary by well...")

df_final['Alignments_Total_SplitAdj'] = \
    df_final.UniqueMappedReads_PreSplit + \
    df_final.UniqueMappedReads_PostSplit / \
    (df_final.NumReadsIn_PostSplit/(df_final.MappingRate_PreSplit*df_final.NumReadsIn))
df_final['MappingRate_Total'] = \
    df_final.Alignments_Total_SplitAdj / df_final.NumReadsIn

# column QC
print("Number of NAs per column:")
print(df_final.isna().sum().to_string())

print("\nNumber of duplicated wells:")
ndupe=df_final.index.duplicated().sum()
print(ndupe)



# final export
df_final.to_csv("Metadata/A07a2_mappingrate.tsv", sep='\t')
print("done.\n\n")
