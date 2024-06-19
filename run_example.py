print('Starting run_example.py')
import os
import pandas as pd
from ALPACA_segment_solution_class import SegmentSolution

# modify this path to run a different cohort, check README for required input files and structure
example_cohort_input_directory = 'input/example_cohort'

model_config = {
    'license': 'remote',
}
config = {'model_config': model_config, 'preprocessing_config': {}}

# if running from bin/MODEL, set wd to ('../../'):
if os.getcwd().endswith('bin/MODEL'):
    print('set working directory to project directory')
    os.chdir('../../')
project_directory = os.getcwd()
cohort_input = example_cohort_input_directory
cohort = cohort_input.split('/')[-1]
assert cohort != ''
# list all directories in the cohort input directory:
tumour_ids = [x for x in os.listdir(cohort_input) if not x.startswith('.')]
tumour_outputs = []
for tumour_id in tumour_ids:
    print(f'running {tumour_id}')
    # set wd to subdirectory of tumour_id with segments
    os.chdir(f'{project_directory}/{cohort_input}/{tumour_id}/segments')
    # get all input files
    input_files = [x for x in os.listdir() if x.endswith('.csv')]
    solutions = []
    for input_file_name in input_files:
        SS = SegmentSolution(input_file_name,config)
        SS.run_iterations()
        SS.find_optimal_solution()
        SS.get_solution()
        solutions.append(SS.optimal_solution)
        print(f'Segment {input_file_name} done.')
    tumour_output = pd.concat(solutions)
    tumour_output_directory = f'{project_directory}/output/{cohort}/patient_outputs/'
    os.makedirs(tumour_output_directory, exist_ok=True)
    tumour_output.to_csv(f'{tumour_output_directory}/final_{tumour_id}.csv', index=False)
    tumour_outputs.append(tumour_output)
    os.chdir(project_directory)
print(f'Creating combined output')
cohort_output = pd.concat(tumour_outputs)
cohort_output_directory = f'{project_directory}/output/{cohort}/cohort_outputs'
os.makedirs(cohort_output_directory, exist_ok=True)
cohort_output.to_csv(f'{cohort_output_directory}/combined.csv', index=False)
print('Done')
print("""
 _____ __    _____ _____ _____ _____
|  _  |  |  |  _  |  _  |     |  _  |
|     |  |__|   __|     |   --|     |
|__|__|_____|__|  |__|__|_____|__|__|
  /\\⌒⌒⌒/\\
  (⦿   ⦿)
  ( 'Y' )
   (   )
   (   )
   (   )
   (~ ~~~~~~~~~~)
   ( ~ ~~   ~~  )
   ( ~  ~ ~  ~  )
   (~  ~~~~~   ~)
    │ │     │ │
    │ │     │ │
""")