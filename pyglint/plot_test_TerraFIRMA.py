import os
import argparse
import json
import yaml
import pickle
import xarray as xr
import matplotlib.pyplot as plt

OUTPUT_DIR = "test_outputs"

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def main():

    ## Get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--icesheet", required=True,
                        help="Ice sheet to plot (e.g. AIS or GrIS)")
    parser.add_argument("--variable", required=True,
                        help="Variable name to plot, must exist in results.nc (e.g. total_smb)")
    args = parser.parse_args()

    icesheet = args.icesheet
    variable_to_plot = args.variable


    ## Load the data from Glint/BISICLES for comparison
    if icesheet == "AIS":
        with open('/gws/ssde/j25b/terrafirma/tm17544/TerraFIRMA_overshoots/processed_data/AIS_data_overview.pkl', 'rb') as file:
            icesheet_d = pickle.load(file) 

        # only keep the first 50 years of data for cs568 to match the pyglint processed data (which only has 50 years of data)
        icesheet_d["cs568"][0] = icesheet_d["cs568"][0].iloc[:10]

    elif icesheet == "GrIS":
        with open('/gws/ssde/j25b/terrafirma/tm17544/TerraFIRMA_overshoots/processed_data/GrIS_data_overview.pkl', 'rb') as file:
            icesheet_d = pickle.load(file)

        # Drop the first year of data as its from an errant file
        icesheet_d["cs568"][0] = icesheet_d["cs568"][0].iloc[1:11]

    # Get the data in the right units/variable for plotting
    bike_total_smb = (icesheet_d["cs568"][0]["grounded_SMB"] + icesheet_d["cs568"][0]["floating_SMB"])*918*1e-12

    plt.figure(figsize=(10, 6))

    for run_folder in sorted(os.listdir(OUTPUT_DIR)):

        run_path = os.path.join(OUTPUT_DIR, run_folder)

        if not os.path.isdir(run_path):
            continue

        config_path = os.path.join(run_path, "config_used.yaml")
        metadata_path = os.path.join(run_path, "metadata.json")
        results_path = os.path.join(run_path, "pyglint_smb_data.nc")

        if not os.path.exists(results_path):
            continue

        cfg = load_yaml(config_path)
        meta = load_json(metadata_path)

        # Filter by icesheet
        if cfg.get("icesheet") != icesheet:
            continue

        # Don't plot masked data for now
        if cfg.get("gr_fl_mask"):
            continue

        # Don't plot cx209 mas data
        if "cx209" in cfg.get("nc_file", ""):
            continue

        # Load timeseries
        ds = xr.load_dataset(results_path)

        if variable_to_plot not in ds:
            continue

        time = ds["time"].values
        values = ds[variable_to_plot].values

        # Build label from metadata + config
        label = (
            f"{meta['git_branch']}"
            #f"{cfg.get('nc_file')[:7]}"
        )

        plt.plot(time[1:], values[1:], label=label)

    # Add Glint/BISICLES data for comparison
    plt.plot(icesheet_d["cs568"][0]["time"]-1850, bike_total_smb, label="Glint/BISICLES", color='black')

    plt.xlabel("Time")
    plt.ylabel(f"{variable_to_plot} (Gt/a)")
    plt.title(f"{icesheet} {variable_to_plot} Comparison")
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig(f"{OUTPUT_DIR}/pyglint_test_TerraFIRMA_{cfg.get('suite_id')}_{variable_to_plot}_{icesheet}.png", dpi=600)

if __name__ == "__main__":
    
    main()