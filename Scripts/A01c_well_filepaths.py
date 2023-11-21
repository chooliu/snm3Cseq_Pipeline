
# ==============================================================================
# Scripts/A01c_well_filepaths.py
# expands plate-level metadata (A01b) into well-level metadata
# ==============================================================================

# recommend running interactively in python/Jupyter to check outputs,
# but shouldn't require any changes to defaults

# load packages ----------------------------------------------------------------

import itertools
import pandas as pd
import numpy as np
import os

# only one potential parameter to change
# by default, processes one row of 384-well plate per "batch"
wells_per_batch = 24 # <--


# # if running interactively, check snm3C_parameters.env loaded or manually spec os.environ e.g.,
# os.environ['dir_proj'] ="/u/home/c/chliu/scratch-cluo/IGVF_iPSC_snm3Cseq_YZCL47"
# os.chdir(os.environ['dir_proj'])
# os.environ['metadat_plate'] = "Metadata/A01b_plate_metadata.csv"
# os.environ['metadat_well'] = "Metadata/A01c_well_filepath.csv"



# expand A01b metadata by well -------------------------------------------------

# load A01b
plates_df = pd.read_csv(os.environ['metadat_plate'], index_col=0)

# from pandas documentation
def expand_grid(data_dict):
    """Create a dataframe from every combination of given values."""
    rows = itertools.product(*data_dict.values())
    return pd.DataFrame.from_records(rows, columns=data_dict.keys())

filepath_df = expand_grid({'plate': plates_df['plate'],
    'row' : [chr(x) for x in range(65, 65+16)],
    'col' : [str(x + 1) for x in range(24)]})
filepath_df['well'] = filepath_df[['row', 'col']].agg(''.join, axis = 1)
filepath_df['wellprefix'] = filepath_df['plate'] + "_" + filepath_df['well']

filepath_df = pd.merge(filepath_df, plates_df, how = "left", on = "plate")



# batch into sets of 24 for bismark mapping, contact calling, etc --------------
# (by default, one row at a time, incremented by "batchnum")

# - alternatively, could make smaller batches of wells (e.g., n = 5) for compute
#   environments that favor many small jobs versus a few long jobs,
# - or... two sets of batches e.g., filepath_df['batchnum_A04a_bismark']
#   pulled by the sub scripts for the A04a script / really resource-intensive jobs only

nwellstot = filepath_df.shape[0]
filepath_df['batchnum'] =\
    pd.Series(range(0, np.ceil(nwellstot / wells_per_batch).astype(int))
             ).repeat(wells_per_batch)[0:nwellstot].reset_index(drop = True) + 1

print( "number of total wells:" )
print( nwellstot )

print( "wells per 'batchnum':" )
print( wells_per_batch )

filepath_df.index = filepath_df.index.astype(int) + 1

def basename(pathin):
    return(pathin.split("/")[-1])

print( "number of plates:" )
print( "Nplates: " + str( filepath_df['platenum'].max() ) )

print( "number of batches:" )
print( "Nbatches: " + str( filepath_df['batchnum'].max() ) )



# then extensive file paths for sections A02-A06 -------------------------------
# (inelegant, but useful for file checking/compiling info)

# A02: demultiplexing 
# all in dir: fastq_demultip/

filepath_df['A02a_fqgz_demultip_R1'] = "fastq_demultip/" + filepath_df[['plate', 'well']].agg('_'.join, axis = 1) + "_indexed_R1.fastq.gz"
filepath_df['A02a_fqgz_demultip_R2'] = "fastq_demultip/" + filepath_df[['plate', 'well']].agg('_'.join, axis = 1) + "_indexed_R2.fastq.gz"

filepath_df['A02a_txt_summary1'] = "fastq_demultip/" + filepath_df['plate'] + "_summary_1.txt"
filepath_df['A02a_txt_summary2'] = "fastq_demultip/" + filepath_df['plate'] + "_summary_2.txt"



# A03: trimming ----------------------------------------------------------------
# all in dir: fastq_trimmed/

filepath_df['A03a_fqgz_paired_R1'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_paired_R1.fastq.gz"
filepath_df['A03a_fqgz_paired_R2'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_paired_R2.fastq.gz"

filepath_df['A03a_fqgz_singletrim_R1'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_singletrim_R1.fastq.gz"
filepath_df['A03a_fqgz_singletrim_R2'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_singletrim_R2.fastq.gz"

filepath_df['A03a_json_fastp'] = "fastq_trimmed/" + filepath_df['wellprefix'] + ".json"


# A04: bismark -----------------------------------------------------------------

filepath_df['A04a_dir_bismark'] = "mapping_bismark/" + filepath_df['wellprefix'] + "/"

# (i) taurus step 1 mapping outputs
filepath_df['A04a_bam_R1p'] = \
    filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_paired_R1'].apply(basename).str.replace(".fastq.gz", "_bismark.bam")
filepath_df['A04a_bam_R2p'] = \
    filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_paired_R2'].apply(basename).str.replace(".fastq.gz", "_bismark.bam")
filepath_df['A04a_bam_R1trims'] = \
        filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_singletrim_R1'].apply(basename).str.replace(".fastq.gz", "_bismark.bam")
filepath_df['A04a_bam_R2trims'] = \
    filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_singletrim_R2'].apply(basename).str.replace(".fastq.gz", "_bismark.bam")

# step 1 logs
filepath_df['A04a_bismarktxt_R1p'] = \
    filepath_df['A04a_dir_bismark'] + filepath_df['wellprefix'] + "_paired_R1_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R2p'] = \
    filepath_df['A04a_dir_bismark'] + filepath_df['wellprefix'] + "_paired_R2_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R1trims'] = \
        filepath_df['A04a_dir_bismark'] + filepath_df['wellprefix'] + "_singletrim_R1_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R2trims'] = \
    filepath_df['A04a_dir_bismark'] + filepath_df['wellprefix'] + "_singletrim_R2_bismark_SE_report.txt"

# (ii) taurus step 2 logs
filepath_df['A04a_bismarktxt_R1p1'] = \
    filepath_df['A04a_dir_bismark'] + "subseq_R1_1_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R1p2'] = \
    filepath_df['A04a_dir_bismark'] + "subseq_R1_2_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R1p3'] = \
    filepath_df['A04a_dir_bismark'] + "subseq_R1_3_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R2p1'] = \
    filepath_df['A04a_dir_bismark'] + "subseq_R2_1_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R2p2'] = \
    filepath_df['A04a_dir_bismark'] + "subseq_R2_2_bismark_SE_report.txt"
filepath_df['A04a_bismarktxt_R2p3'] = \
    filepath_df['A04a_dir_bismark'] + "subseq_R2_3_bismark_SE_report.txt"

# (iii) picard de-duplication
filepath_df['A04a_log_picard'] = filepath_df['A04a_dir_bismark'] + "picard.log"

# (iv) final merged, sorted, dedupe bam
filepath_df['A04a_bam_final'] = filepath_df['A04a_dir_bismark'] + "merged_dedupe.bam"

# A04c: mapping stats ----------------------------------------------------------

filepath_df['A04c_txt_samstats'] = filepath_df['A04a_dir_bismark'] + "samstats.txt"
filepath_df['A04c_txt_covtot'] = filepath_df['A04a_dir_bismark'] + "nbases_cov_by_chr.txt"
filepath_df['A04c_txt_covnsites'] = filepath_df['A04a_dir_bismark'] + "total_cov_by_chr.txt"

# A05a: methylation quantification ---------------------------------------------

filepath_df['A05a_allc'] = filepath_df['A04a_dir_bismark'] + "allc.tsv.gz"
filepath_df['A05a_allctbi'] = filepath_df['A04a_dir_bismark'] + "allc.tsv.gz.tbi"

# A06: contact mapping ---------------------------------------------------------

filepath_df['A06a_pairs'] = filepath_df['A04a_dir_bismark'] + "pairs.tsv"
filepath_df['A06a_3c_metadat'] = filepath_df['A04a_dir_bismark'] + "metadat_pairs.tsv"



# finally, export --------------------------------------------------------------
# by default exports to Metadata/A01c_well_filepath.csv


print("final metadata file dimensions:")
print(filepath_df.shape)
filepath_df.to_csv(os.environ['metadat_well'])



