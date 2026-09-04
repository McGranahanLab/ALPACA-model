import numpy as np
import pandas as pd


def get_chr_table(chr_table_path):
    chr_table = pd.read_csv(chr_table_path)
    chr_table['cumsum'] = np.cumsum(chr_table['len'])
    chr_table['shift'] = [0] + list(chr_table['cumsum'][:-1])
    chr_table['ticks'] = chr_table['shift'] + chr_table['len'] / 2
    return chr_table


def remove_duplicates_preserve_order(input_list):
    seen = {}
    result = []
    for item in input_list:
        if item not in seen:
            seen[item] = True
            result.append(item)
    return result


def get_tree_edges(tree_paths):
    all_edges = list()
    for path in tree_paths:
        if len(path) == 2:
            all_edges.append(tuple(path))
        else:
            for i in range(len(path) - 1):
                all_edges.append((path[i], path[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges


def get_unique_lists(list_of_lists):
    u = []
    for e in list_of_lists:
        if e not in u:
            u.append(e)
    return u


def clean_segment(seg):
    if isinstance(seg, str):
        if seg.isdigit():
            seg = int(seg)
        else:
            seg = np.nan
    return seg


def clean_output(patient_output):
    if ('abs_segment' in patient_output.columns) and ('segment' not in patient_output.columns):
        patient_output = patient_output.T.drop_duplicates().T
        patient_output['segment'] = patient_output['abs_segment']
    if len(patient_output["segment"].iloc[0].split('_')) == 3:
        patient_output["chr"] = patient_output["segment"].str.split("_", expand=True)[0].str.replace('chr', '').apply(lambda x: clean_segment(x))
        patient_output['start'] = patient_output["segment"].str.split("_", expand=True)[1].apply(lambda x: clean_segment(x))
        patient_output['end'] = patient_output["segment"].str.split("_", expand=True)[2].apply(lambda x: clean_segment(x))
    else:  # this means that chr column is already unpacked and that segment is in format 123_345 (no chr)
        patient_output['start'] = patient_output["segment"].str.split("_", expand=True)[0].apply(lambda x: clean_segment(x))
        patient_output['end'] = patient_output["segment"].str.split("_", expand=True)[1].apply(lambda x: clean_segment(x))
    patient_output = patient_output.loc[~(patient_output["chr"].isna() | (patient_output["start"].isna()) | (patient_output["end"].isna()))]
    return patient_output


def find_parent(clone, tree):
    for branch in tree:
        if clone in branch:
            clone_index = branch.index(clone)
            if clone_index != 0:
                return branch[clone_index - 1]
            else:
                return 'diploid'
