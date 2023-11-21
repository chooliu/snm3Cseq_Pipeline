
# A06a_quantify_contacts_TAURUS.py, v0.1 =======================================
# this is a re-code of TAURUS-MH by @chooliu 
# trades off some readability/modularity versus efficiency (reads in memory)
# ==============================================================================

import sys
import os
import pandas as pd
import numpy as np
import bamread



# define valid chromosomes =====================================================

# here just for for chromosome ordering / order of pairs in final output
# select chromosomes to include in contact matrices

chrom_sizes = \
    pd.read_csv(os.environ['ref_chromsizes'], sep="\t", header=None
               ).set_axis(['chr', 'len'], axis = 1)
chrom_sizes['chr'] = pd.Categorical(chrom_sizes['chr'],
    categories = chrom_sizes['chr'] , ordered = True)

# by default include anything in genome ref files...
list_valid_chrom = chrom_sizes['chr']

# # but some example alternatives below: could be all autosomal contacts, or +XYM
# # caution: check how yr downstream HiC/3C software of choice deals with non-autosomal chr
# list_valid_chrom = ["chr" + str(i) for i in range(1, 99)] + ["chrX", "chrY", "chrM"]
# list_valid_chrom = ["chr" + str(i) for i in range(1, 99)] 

def tidy_chr_order(x):
    return(pd.Categorical(x, categories = list_valid_chrom, ordered = True))



# load alignments & process info per read name =================================
# check format, fields, filtering (e.g., add MAPQ filt) if aligner changes

alignments = bamread.read_bam_full(sys.argv[1])

alignments = \
    alignments.drop(["QueryStart", "Cigar", "Quality", "QuerySequence"], axis = 1
               ).rename({'QueryEnd': 'Length'}, axis = 1)

# parse .fastq file names (including edits from 'split')
# order of 'split' category represents preference when taking "outer-most" contacts
alignments['ReadPrefix'] = alignments['Name'].str.split(":N:0:").str[0]
alignments['Split'] = alignments['ReadPrefix'].str.split("_").str[1] 
alignments['Split'] = pd.Categorical(
    alignments['Split'],
    categories=["1", "1:P1", "1:P2", "1:P3", "2:P3", "2:P2", "2:P1", "2"],
    ordered=True)
alignments['ReadPrefix'] = alignments['ReadPrefix'].str.split("_").str[0]
alignments['ReadPair'] = alignments['Split'].str.split(":").str[0] # 1 or 2

# 'pos' represents the 5' position of each read
# depending on convention could also report "leftmost" (smallest val on REF) or "rightmost"
alignments['Pos'] = alignments['Start'].where(
    alignments['ReadPair']=="1", other=alignments['End'])

# let a candidate contact = read names with more than
# read/split alignment type -- (1, 2, 1:P1, ..., 2:P3)
readprefix_candidate_contact = alignments.groupby("ReadPrefix")['Split'].nunique()

# distribution of splits aligned per readname
# possible values are between 1 and 4
n_align_per_read = readprefix_candidate_contact.value_counts().sort_index()

# examine readprefixes with n >= 2 splits per readname
# as they are candidates that contain possible intra- or intra-chromosomal info about contacts
readprefix_candidate_contact = readprefix_candidate_contact[readprefix_candidate_contact >= 2].index

# define potential contacts as "outermost" read 5'-positions a la TAURUS-MH
# among readprefixes with multiple alignments
# (only one contact possible per read pair)
contacts_taurus = \
    alignments[alignments['ReadPrefix'].isin(readprefix_candidate_contact)
              ].groupby(['ReadPrefix']
              ).filter(lambda x : len(x) > 1
              ).groupby(['ReadPrefix']).nth([0, -1]
              ).reset_index()

# for readprefixes with intrachr, sorts so that chr1 < chr3 first in df
contacts_taurus = \
    contacts_taurus.assign(Chromosome=tidy_chr_order(contacts_taurus.loc[:, 'Chromosome'])
                          ).sort_values(['ReadPrefix', 'Chromosome', 'Pos'])


# note: may keep read name in future (to be true .pairs format)
# would need to add new column with readname
contacts_taurus = contacts_taurus.drop('ReadPrefix', axis = 1)



# TAURUS-like --> contacts file ================================================

# initialize
pairs_taurus = \
    pd.DataFrame(index = range(int((contacts_taurus.shape[0] + 1)/2)),
                 columns = range(6))
pairs_taurus.columns = ['chr1', 'pos1', 'chr2', 'pos2', 'strand1', 'strand2']

# alignments --> pairs
pairs_taurus.iloc[:, [0, 1]] = \
    contacts_taurus[['Chromosome', 'Pos']].iloc[0::2, ].reset_index(drop = True)
pairs_taurus.iloc[:, [2, 3]] = \
    contacts_taurus[['Chromosome', 'Pos']].iloc[1::2, ].reset_index(drop = True)
pairs_taurus.iloc[:, 4] = contacts_taurus['Strand'].iloc[0::2, ]
pairs_taurus.iloc[:, 5] = contacts_taurus['Strand'].iloc[1::2, ]

# for intra distances, calculate distance
pair_is_intra = pairs_taurus['chr1'] == pairs_taurus['chr2']
intra_dist = pairs_taurus['pos2'] - pairs_taurus['pos1']

# optional: exclude close by contacts that may be re-ligations or non-contact (e.g., >1kb only)
# note: there can small # pairs with intra dist exactly == 0 excluded here (short insert sizes)
min_intra_dist = 0 # <-- optionally increase from default of 0
filt_intra_dist = intra_dist > min_intra_dist

# duplicated chr pos strand (both pairs)
filt_nondupe = ~pairs_taurus.duplicated()

# chromosome valid
filt_chrom = pairs_taurus['chr1'].isin(list_valid_chrom) & \
    pairs_taurus['chr2'].isin(list_valid_chrom)



# print out final =======================================================================

# criteria for final filtering:
# - either intra-chrom greater than minimum distance threshold, or inter-chrom
# - within target chromosomes ("list_valid_chrom")
# - non-duplicated by position
filter_final = (~ pair_is_intra | filt_intra_dist ) & filt_chrom & filt_nondupe
final_out = pairs_taurus[filter_final]

final_out = final_out.assign(chr1=tidy_chr_order(final_out.loc[:, 'chr1'])
                            ).assign(chr2=tidy_chr_order(final_out.loc[:, 'chr2']))
final_out = final_out.sort_values(['chr1', 'chr2', 'pos1', 'pos2'])

final_out.to_csv("pairs.tsv", sep = "\t", index_label=False, index = False, header=False)



# final output stats =======================================================================

dict_output_stats = {

    # readnames
    'readnames_anyalignment' : alignments['ReadPrefix'].nunique(),     # readnames with at least one alignment, in any split type (1, 1:P1, ..., 2:P3, 2)
    'readnames_singlealign' : n_align_per_read.loc[1, ], # read names with only one alignment (non-informative)
    
    # alignments
    'alignments_total' : alignments.shape[0],     # total number of alignments
    'alignments_wholeread' : sum(alignments['Split'].isin(['1', '2'])),     # alignments based on 1 or 2
    'alignments_taurus_3split' : sum(~alignments['Split'].isin(['1', '2'])), # alignments from 1:P1, ..., 2:P3

    # contacts prefiltering
    'contacts_prefilt' : len(readprefix_candidate_contact), # num candidate contacts pre-filt
    'contacts_prefilt_minintralen' : sum(~(~ pair_is_intra | filt_intra_dist )), # inter or >min_intra_dist
    'contacts_prefilt_nondupe' : sum(~filt_nondupe), # duplicated chr:pos:strand
    'contacts_prefilt_validchrom' : sum(~filt_chrom), # "list_valid_chrom"
    
    # optional pre-filtering info --
    # comment in below if want to check any bias induced by "pre-filtering"
#     'contacts_prefilt_intra' : sum(pair_is_intra), # pair is on same chrom
#     'contacts_prefilt_inter' : sum(~pair_is_intra), # pair on different chromosomes
#     'contacts_prefilt_intra_less1000' : sum(pair_is_intra & (intra_dist <= 1000)), # same chromosome
#     'contacts_prefilt_intra_greater1000' : sum(pair_is_intra & (intra_dist > 1000)), # same chrom, >1k
#     'contacts_prefilt_intra_greater10000' : sum(pair_is_intra & (intra_dist > 10000)), # same chrom >10k
    
    # final contacts
    'contacts_final' : sum(filter_final), # final number contacts
    'contacts_final_intra' : sum(filter_final & pair_is_intra), # final, same chrom
    'contacts_final_inter' : sum(filter_final & ~pair_is_intra), # final, diff chrom

    'contacts_final_intra_less1000' : sum(filter_final & pair_is_intra & (intra_dist <= 1000)), # final, diff chrom
    'contacts_final_intra_greater1000' : sum(filter_final & pair_is_intra & (intra_dist > 1000)), # final, diff chrom
    'contacts_final_intra_greater10000' : sum(filter_final & pair_is_intra & (intra_dist > 10000))

}



# fractions / ratios of values calculated above

dict_output_stats['frac_readnames_yieldcontact'] = \
    dict_output_stats['contacts_final']/dict_output_stats['readnames_anyalignment']
dict_output_stats['frac_contacts_passfilt'] =\
    dict_output_stats['contacts_final']/dict_output_stats['contacts_prefilt']
    
dict_output_stats['frac_finalcont_inter'] = \
    dict_output_stats['contacts_final_inter']/dict_output_stats['contacts_final']
dict_output_stats['frac_finalcont_intraless1kb_totintra'] = \
    dict_output_stats['contacts_final_intra_less1000']/dict_output_stats['contacts_final_intra']
dict_output_stats['frac_finalcont_intragreater1kb_totintra'] = \
    dict_output_stats['contacts_final_intra_greater1000']/dict_output_stats['contacts_final_intra']
dict_output_stats['frac_finalcont_intragreater10kb_totintra'] = \
    dict_output_stats['contacts_final_intra_greater10000']/dict_output_stats['contacts_final_intra']
    
dict_output_stats['ratio_finalcont_intragreater1kb_inter'] = \
    dict_output_stats['contacts_final_intra_greater1000']/dict_output_stats['contacts_final_inter']


# final output
pd.DataFrame.from_dict(dict_output_stats, orient = "index"
    ).transpose().to_csv("metadat_pairs.tsv", sep = "\t", index_label=False, index = False)
    
print("finished processing " + sys.argv[1] + ".")
