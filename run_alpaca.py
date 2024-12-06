#!/usr/bin/env python3
import argparse
import shlex
import time
from ALPACA_segment_solution_class import SegmentSolution

parser = argparse.ArgumentParser(description="Run ALPACA with specified parameters.")
parser.add_argument('--input_data_directory',type=str,required=True,help="Directory where input data is stored. Should contain subdirectories for each tumour")
parser.add_argument('--use_binary_search',type=int,default=False,help="Whether to use binary search or not. Binary search is faster but may not find optimal solution.")
parser.add_argument('--use_two_objectives',type=int,default=False,help="Whether to use two objectives or not. First objective minimises number of segments outside CI, second objective minimises error.")
parser.add_argument('--use_minimise_events_to_diploid',type=int,default=False,help="Whether to minimize events to diploid or not. If true, ALPACA will introduce diploid pseudo-clone at the root of the tree.")
parser.add_argument('--exclusive_amp_del',type=int,default=1,help="")
parser.add_argument('--input_files',type=str,required=True,help="Input table for one segment.")
parser.add_argument('--prevent_increase_from_zero_flag', type=int, required=True,help='')
parser.add_argument('--add_event_count_constraints_flag', type=int, required=True,help='')
parser.add_argument('--add_state_change_count_constraints_flag', type=int, required=True,help='')
parser.add_argument('--add_path_variability_penalty_constraints_flag', type=int, required=True,help='')
parser.add_argument('--add_allow_only_one_non_directional_event_flag', type=int, required=True,help='')
parser.add_argument('--time_limit', default=60, type=int, help='Time limit in seconds for each model run')
parser.add_argument('--homozygous_deletion_threshold',type=float, help='')
parser.add_argument('--missing_clones_inherit_from_children_flag', default=1, type=int, help='Ensure that missing clones inherit cn from childen (events go up in the tree)')
parser.add_argument('--cpus', default=2, type=int, help='number of available cpus') 
parser.add_argument('--rsc', default=0, type=int, help='remove small clones') 
parser.add_argument('--ccp', default=0, type=int, help='calibrate clone proportions') 
parser.add_argument('--ci', default=0.9, type=float, help='Confidence interval for SNP copynumbers')
parser.add_argument('--ci_table_name', default='', type=str, help='Name of file containing confidence intervals for SNP copynumbers')
parser.add_argument('--s_type', default='s_strictly_decreasing', type=str)
parser.add_argument('--compare_with_true_solution', default=0, type=int)
parser.add_argument('--run_with_qc', default=0, type=int)
parser.add_argument('--output_all_solutions', default=0, type=int)
parser.add_argument('--output_model_selection_table', default=0, type=int)
# parse arguments
args = parser.parse_args()
input_files = shlex.split(args.input_files)  # each input file is a file name of a table containing single segment

# make config dictionary
model_config = {
    'use_binary_search': args.use_binary_search,
    'use_two_objectives': args.use_two_objectives,
    'use_minimise_events_to_diploid': args.use_minimise_events_to_diploid,
    'exclusive_amp_del': args.exclusive_amp_del,
    'prevent_increase_from_zero_flag': args.prevent_increase_from_zero_flag,
    'add_event_count_constraints_flag': args.add_event_count_constraints_flag,
    'add_state_change_count_constraints_flag': args.add_state_change_count_constraints_flag,
    'add_path_variability_penalty_constraints_flag': args.add_path_variability_penalty_constraints_flag,
    'add_allow_only_one_non_directional_event_flag': args.add_allow_only_one_non_directional_event_flag,
    'homozygous_deletion_threshold':args.homozygous_deletion_threshold,
    'homo_del_size_limit': 5 * 10 ** 7,
    'time_limit': args.time_limit,
    'cpus': args.cpus}
preprocessing_config = {
    'rsc': args.rsc,
    'ccp': args.ccp,
    'ci': args.ci,
    'ci_table_name': args.ci_table_name,
    'compare_with_true_solution': args.compare_with_true_solution,
    's_type': args.s_type,
    'input_data_directory': args.input_data_directory,
    'run_with_qc': args.run_with_qc,
    'output_all_solutions': args.output_all_solutions,
    'output_model_selection_table' : args.output_model_selection_table
}

config = {'model_config': model_config, 'preprocessing_config': preprocessing_config}


print('-------------------------------------------------')
print('Running ALPACA with the following parameters:')
# print value of each parameter:
for arg in vars(args):
    print(arg, getattr(args, arg))


for input_file_name in input_files:
    output_name = 'optimal'+input_file_name.split('ALPACA_input_table_')[1]
    start_time = time.time()
    SS = SegmentSolution(input_file_name, config)
    SS.run_iterations()
    SS.find_optimal_solution()
    SS.get_solution()
    optimal_solution = SS.optimal_solution
    assert optimal_solution is not None
    if SS.output_all_solutions:
        all_solutions = SS.get_all_simplified_solution()
        all_solutions_output_name = output_name.replace('optimal','all')
        all_solutions.to_csv(all_solutions_output_name, index=False) # TODO change to parquet
    if SS.output_model_selection_table:
        output_model_selection_table = SS.elbow_search_df_strictly_decreasing
        output_model_selection_table.to_csv(f'{SS.tumour_id}_{SS.segment}_model_selection_table.csv', index=False)
    end_time = time.time()
    total_run_time = round(end_time - start_time)
    optimal_solution['run_time_seconds'] = total_run_time
    optimal_solution.to_csv(output_name, index=False)
    print(f'Segment {input_file_name} solved.')
