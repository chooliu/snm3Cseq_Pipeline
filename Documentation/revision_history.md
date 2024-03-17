
# Revision History
Key updates from past iterations (& some superfluous details as notes for myself)

# Public Release, v1

### v1.0.1
- Contact pairs now outputted as gzip'd file (`pairs.tsv.gz`)
- Added `Documentation/submission_helper.txt`
- Minor bug fix to creation of contact pairs file in `Scripts/A06a_quantify_contacts_TAURUS.py`, in which the strand column was intermittently missing, although this wasn't used in most downstream applications. (potential `reset_index()` issue, and `pandas` "DeprecationWarning: In a future version, `df.iloc[:, i] = newvals` will attempt to set the values inplace instead of always setting a new array")
- Added convenience parsing of .env file variables for some interactive scripts 

### v1.0.0
- First public release of `snm3Cseq_taurus` based on TAURUS-MH.
- Originally released Nov 20th, 2023; two minor bug fixes (without impact on output) and finished documentation writing on Nov 29th, 2023.
