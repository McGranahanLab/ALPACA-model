import os
import sys
sys.path.append('.')
from ALPACA_segment_solution_class import SegmentSolution
from run_alpaca import solve
import pandas as pd

tumour_id = 'SA1049'
cohort = 'sc_psedubulk_thr_30'
output_dir = '/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/publication/output/sc_psedubulk_thr_30_interactive_test/cohort_outputs'
os.makedirs(output_dir, exist_ok=True)

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
    'ci_table_name':'',
    'output_all_solutions':0
}

config = {'model_config': model_config, 'preprocessing_config': preprocessing_config}


os.chdir(f'../../input/{cohort}/{tumour_id}/segments')

all_segments = [x.replace('.csv','').split(f'{tumour_id}_')[1] for x in os.listdir() if x.endswith('csv')]

solutions = []
for segment in all_segments[0:4]:
    input_file_name = f'ALPACA_input_table_{tumour_id}_{segment}.csv'
    SS = SegmentSolution(input_file_name, config)
    SS.run_iterations()
    SS.find_optimal_solution()
    # SS.compare_with_true_solution()
    SS.get_solution()
    optimal_solution = SS.optimal_solution
    if SS.output_all_solutions:
        all_solutions = SS.get_all_simplified_solution()
    print(f'{segment} done!')
    solutions.append(optimal_solution)
solutions = pd.concat(solutions)
solutions.to_csv(f'{output_dir}/combined.csv', index=False)