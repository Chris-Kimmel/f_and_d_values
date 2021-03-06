#!/usr/bin/env python3
'''
f_and_d_values.py

Compute the fractions (and dampened fractions) of modifed reads from a p-value
CSV file.

Written for the Kim Lab (Aug 9, 2021)

### Usage

Example:
./f_and_d_values.py /fs/project/PAS1405/GabbyLee/project/m6A_modif/WT_cellular/23456_WT_cellular.csv output.csv

This script takes two positional arguments from the command line, which are, in
order:
    - filepath_to_read 
    - filepath_to_write

The user may also wish to change the values of the constants LOWER_THRESH and
UPPER_THRESH, which are hardcoded in this script.

### Other

The results of this script will be close to the values "tombo text_output
browser_files" would produce, but they will still differ slightly.

Chris Kimmel
chris.kimmel@live.com
'''

################################################################################
############################# imports and pragmas ##############################
################################################################################

# pylint: disable=invalid-name,line-too-long

import pandas as pd
from sys import argv


################################################################################
################################## constants ###################################
################################################################################

# These are constants in Tombo's f-value and d-value formulas.
LOWER_THRESH = 0.05 # Tombo's default value for RNA is 0.05 for this variable.
UPPER_THRESH = 0.40 # Tombo's default value for RNA is 0.40 for this variable.


################################################################################
############################ command-line arguments ############################
################################################################################

# full filepaths
filepath_to_read = argv[1]
filepath_to_write = argv[2]


################################################################################
################################# subroutines ##################################
################################################################################

def load_csv(filepath):
    '''Load per-read stats from a CSV file into a Pandas DataFrame'''
    retval = (
        pd.read_csv(filepath, header=0, index_col=0)
        .rename_axis('pos_0b', axis=1)
        .rename_axis('read_id')
    )
    retval.columns = retval.columns.astype(int)
    return retval

def longify(df):
    '''Convert dataframe output of load_csv to a long format'''
    return df.stack().rename('pval').dropna().reset_index()

# widify is not called in this script
def widify(df):
    '''Undo longify'''
    return df.set_index(['read_id', 'pos_0b']).unstack().droplevel(0, axis=1)


################################################################################
################################## load data ###################################
################################################################################

wide = load_csv(filepath_to_read)

if not wide.index.is_unique:
    raise NotImplementedError("There is a duplicate read ID in the given file. "
        "This script requires unique read IDs.")

raw = longify(wide)
del wide


################################################################################
############################ compute d and f values ############################
################################################################################

results = (
    raw
    .assign(below_lower_thresh=lambda x: x['pval'] < LOWER_THRESH)
    .assign(above_upper_thresh=lambda x: x['pval'] > UPPER_THRESH)
    .groupby('pos_0b')
    .agg(
        num_below_lower_thresh=('below_lower_thresh', sum),
        num_above_upper_thresh=('above_upper_thresh', sum),
        covg=('read_id', 'count'),
    )
    .reset_index()
    .assign(
        frac_below_lower_thresh=lambda x: x['num_below_lower_thresh'] / x['covg'],
        frac_above_upper_thresh=lambda x: x['num_above_upper_thresh'] / x['covg'],
    )
    .assign(
        f_value=lambda x: x['num_below_lower_thresh'] / (
            x['num_below_lower_thresh'] + x['num_above_upper_thresh'] ),
        d_value=lambda x: x['num_below_lower_thresh'] / (
            x['num_below_lower_thresh'] + x['num_above_upper_thresh'] + 2 ),
    )
)


################################################################################
################################## write data ##################################
################################################################################

results.to_csv(filepath_to_write, index=False)
