#!/usr/bin/env python3
import sys
from pathlib import Path

try:
    from importlib.metadata import version, PackageNotFoundError
except Exception:
    try:
        from importlib_metadata import version, PackageNotFoundError  # type: ignore
    except Exception:
        # If importlib.metadata isn't available, provide a fallback function
        def version(pkg_name):
            return "unknown"


def _print_version():
    try:
        v = version("alpaca")
    except Exception:
        v = "unknown"
    print(v)


def main():
    command = sys.argv[1].lower() if len(sys.argv) > 1 else ""

    if command in ["version", "--version", "-v"]:
        _print_version()
        return

    if command in ["help", "--help", "-h"]:
        from alpaca.utils import show_help

        show_help()
        return

    if command == "input-conversion":
        import alpaca.scripts as scripts

        rc = scripts.input_conversion()
        if rc:
            sys.exit(rc)
        return
    elif command == "ancestor-delta":
        import alpaca.scripts as scripts

        scripts.run_get_cn_change_to_ancestor()
        return
    elif command == "ccd":
        import alpaca.scripts as scripts

        scripts.run_calculate_ccd()
        return
    elif command == "wgd":
        import alpaca.scripts as scripts

        scripts.run_calculate_wgd()
        return
    elif command == "plot-tumour":
        import alpaca.scripts as scripts

        scripts.run_plot_tumour()
        return
    elif command == "run":
        run_alpaca()
        return
    elif command == "get-scores":
        import alpaca.scripts as scripts

        scripts.run_get_scores()
        return

    print(f"Unknown command: {command}")
    print("Run 'alpaca help' for available commands.")
    sys.exit(1)


def run_alpaca():
    import os
    from tqdm import tqdm
    from io import StringIO

    from alpaca.ALPACA_segment_solution_class import SegmentSolution
    from alpaca.utils import (
        print_logo,
        concatenate_output,
        set_run_mode,
        create_logger,
        save_dataframe_to_csv,
        process_ci_reports,
        process_monoclonal_reports,
        process_run_summary_reports,
        process_infeasibility_reports,
        SegmentInfeasibleError,
    )
    from alpaca.make_configuration import make_config
    from alpaca.analysis import get_cn_change_to_ancestor

    # Configure logging
    logger = create_logger(name="ALPACA", log_dir="logs")
    logger.info("Starting ALPACA")
    config = make_config(sys.argv[1:])
    debug = config["preprocessing_config"]["debug"]
    output_all = config["preprocessing_config"].get("output_all_solutions", False)
    if output_all:
        out_dir = config["preprocessing_config"].get("output_directory", "./")
        all_solutions_dir = os.path.join(out_dir, "all_solutions")
        try:
            os.makedirs(all_solutions_dir, exist_ok=True)
            logger.info(f"Created directory for all solutions: {all_solutions_dir}")
        except Exception as e:
            logger.error(f"Failed to create all_solutions directory {all_solutions_dir}: {e}")
            raise
    if debug:
        logger.setLevel("DEBUG")
        logger.info("Debug mode is ON")
        # to enable testing model with a provided solution, check if the solution dataframe is provided.
        # we expect only single segmnt there, so throw an error if there are more.
        # then, modify 'input_files' to contain only this single segment file.
        if config["preprocessing_config"].get("test_with_provided_solution") is not None:
            import pandas as pd
            provided_solution = pd.read_csv(
                config["preprocessing_config"]["test_with_provided_solution"]
            )
            segments = provided_solution["segment"].unique()
            if len(segments) > 1:
                raise ValueError(
                    "When using 'test_with_provided_solution', the provided solution file must contain only a single segment."
                )
            config["preprocessing_config"]["input_files"] = [x for x in config["preprocessing_config"]["input_files"] if segments[0] in os.path.basename(x)]
            
    # determine running mode:
    # if 'tumour', expect single file with all the segments and output a single file
    # if 'segment' expect array of files to segment files (can be from different tumours) and create separate outputs for each segment
    config, run_mode = set_run_mode(config)
    # if we run in 'segment' mode, disable plotting by default:
    if run_mode == "segment":
        config["preprocessing_config"]["plot_output_mode"] = "none"
    logger.info("-------------------------------------------------")
    logger.info("Running ALPACA with the following parameters:")
    # print value of each parameter:
    logger.info(config)
    print_logo()
    # initiate progress bar:
    if not debug:
        progress_bar = tqdm(
            total=len(config["preprocessing_config"]["input_files"]),
            desc="Processing files",
            unit="file",
            file=sys.stderr,
        )
        original_stdout = sys.stdout
        if os.name == "nt":  # Windows
            sys.stdout = open(os.devnull, "w")
        else:  # Unix/Linux
            sys.stdout = StringIO()
    try:
        for input_file_name in config["preprocessing_config"]["input_files"]:
            SS = SegmentSolution(input_file_name, config, logger)
            if (
                not config["preprocessing_config"]["overwrite_output"]
                and SS.output_exists()
            ):
                logger.warning(
                    f"Output for {input_file_name} already exists. Use '--overwrite_output 1' option to overwrite existing output. Skipping this segment."
                )
                continue
            logger.debug(f"Output path: {SS.create_output_path()}")
            try:
                SS.run_iterations()
                SS.find_optimal_solution()
                SS.get_solution()
                SS.save_output()
            except SegmentInfeasibleError as e:
                output_dir = config["preprocessing_config"].get("output_directory", ".")
                logger.warning(
                    f"Segment infeasible — skipping. tumour={e.tumour_id}, "
                    f"segment={e.segment}, solver_status={e.solver_status}. "
                    "An infeasibility report has been written to the output directory."
                )
                if not debug:
                    progress_bar.update(1)
                    progress_bar.set_description(f"INFEASIBLE {input_file_name}")
                continue
            if not debug:
                progress_bar.update(1)
                progress_bar.set_description(f"Processing {input_file_name}")
            else:
                logger.info(f"Segment {input_file_name} solved.")
        if run_mode == "tumour":
            concatenated_output_path = concatenate_output(
                config["preprocessing_config"]["output_directory"]
            )
            logger.info("Calculating copy number change to ancestor...")
            cn_change_to_ancestor_df = get_cn_change_to_ancestor(
                f"{SS.tumour_dir}/tree_paths.json", concatenated_output_path
            )
            save_dataframe_to_csv(
                df=cn_change_to_ancestor_df,
                output_dir=SS.config["preprocessing_config"]["output_directory"],
                output_filename="cn_change_to_ancestor.csv",
            )
            # parse and combine reports:
            output_dir = SS.config["preprocessing_config"]["output_directory"]
            process_ci_reports(output_dir, delete=True, outpath=output_dir + "/ci_modified_report.csv")
            process_monoclonal_reports(output_dir, delete=True, outpath=output_dir + "/monoclonal_samples_report.csv")
            process_run_summary_reports(output_dir, delete=True, outpath=output_dir + "/run_gap_summary.csv")
            process_infeasibility_reports(output_dir, delete=True, outpath=output_dir + "/infeasibility_report.csv")
            # cleanup: move all the logs, reports and summaries to a separate folder:
            logs_and_reports_dir = os.path.join(output_dir, "logs_and_reports")
            try:
                import shutil
                os.makedirs(logs_and_reports_dir, exist_ok=True)
                for item in os.listdir(output_dir):
                    if item == "logs_and_reports":
                        continue
                    if "_log_" in item or "report" in item or item.endswith("_summary.csv"):
                        shutil.move(os.path.join(output_dir, item), os.path.join(logs_and_reports_dir, item))
                logger.info(f"Moved logs and reports to {logs_and_reports_dir}")
            except Exception as e:
                logger.error(f"Failed to move logs and reports to {logs_and_reports_dir}: {e}")
                raise
            plot_mode = config["preprocessing_config"].get("plot_output_mode", "notebook")
            heatmap_palette = config["preprocessing_config"].get("heatmap_palette")
            genome_build = config["preprocessing_config"].get("genome_build", "hg19")
            if plot_mode != "none":
                from alpaca.plotting import export_plot_outputs

                tumour_input_dir = config["preprocessing_config"].get("input_tumour_directory") or SS.tumour_dir
                try:
                    export_plot_outputs(
                        mode=plot_mode,
                        tumour_input_dir=Path(tumour_input_dir),
                        tumour_output_dir=Path(output_dir),
                        alpaca_output_path=Path(concatenated_output_path),
                        heatmap_palette=heatmap_palette,
                        genome_build=genome_build,
                    )
                    logger.info(
                        f"Generated {plot_mode} plot artefacts in {output_dir}"
                    )
                except Exception as plot_error:
                    logger.error(f"Failed to generate plots ({plot_mode}): {plot_error}")
                    if debug:
                        raise
            logger.info(
                f"""Analysis completed successfully. Output saved to: {output_dir}"""
            )
        logger.info("Done")
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise e
    finally:
        if not debug:
            sys.stdout.close()
            sys.stdout = original_stdout
            progress_bar.close()


if __name__ == "__main__":
    main()
