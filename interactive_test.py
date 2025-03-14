import os
import sys
sys.path.append('.')
from ALPACA_segment_solution_class import SegmentSolution

from alpaca.run_alpaca import solve
tumour_id = 'LTXSIM127'
segment = '10_1018473_1437135'
cohort = 'simulations'
input_file_name = f'ALPACA_input_table_{tumour_id}_{segment}.csv'


model_config = {
    'use_binary_search': 0,
    'use_two_objectives': 1,
    'use_minimise_events_to_diploid': 1,
    'exclusive_amp_del': 1,
    'prevent_increase_from_zero_flag': 1,
    'add_event_count_constraints_flag': 1,
    'add_state_change_count_constraints_flag': 0,
    'add_path_variability_penalty_constraints_flag': 0,
    'add_allow_only_one_non_directional_event_flag': 1,
    'homozygous_deletion_threshold': 1,
    'homo_del_size_limit': 5 * 10 ** 7,
    'variability_penalty': 1,
    'time_limit': 60,
    'cpus': 2
}

preprocessing_config = {
    'rsc': False,
    'ccp': True,
    'ci': 0.5,
    'compare_with_true_solution': False,
    's_type': 's_strictly_decreasing',
    'chr_table_file': '../_assets/chr_table.csv',
    'input_data_directory': f'../input/{cohort}',
    'ci_table_name':'ci_table.csv',
    'run_with_qc':0,
    'output_all_solutions':0,
    'output_model_selection_table':1
}

config = {'model_config': model_config, 'preprocessing_config': preprocessing_config}


os.chdir(f'../../input/{cohort}/{tumour_id}/segments')


SS = SegmentSolution(input_file_name, config)
SS.run_iterations()
SS.find_optimal_solution()
# SS.compare_with_true_solution()
SS.get_solution()
optimal_solution = SS.optimal_solution
if SS.output_all_solutions:
    all_solutions = SS.get_all_simplified_solution()
all_solutions = SS.get_all_simplified_solution()
elbow_table = SS.elbow_search_df_strictly_decreasing
print('done!')