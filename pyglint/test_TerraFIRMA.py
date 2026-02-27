# Imports for config management
import argparse
from datetime import datetime
import os
import subprocess
import sys
import json
import yaml

# Imports for testing pyglint
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
import xarray as xr

from coupling import atm_to_ism
from transformation import up_down_pair, Uniform2DGrid
from transformation import PROJ_ARCTIC_4326, PROJ_ANTARCTIC_3031

# Utility: load YAML config
def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)
    

# Utility: get Git info
def get_git_info():
    try:
        commit = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD']
        ).decode().strip()
        branch = subprocess.check_output(
            ['git', 'rev-parse', '--abbrev-ref', 'HEAD']
        ).decode().strip()
        return {'git_commit': commit, 'git_branch': branch}
    except Exception:
        return {'git_commit': 'unknown', 'git_branch': 'unknown'}
    

# Get ice sheet grid and mask
def read_ism_nc(cfg):

    nc_file = cfg['nc_file']
    nc_bike_in = Dataset(nc_file,'r')
    
    if cfg['icesheet'] == 'GrIS':

        topo = nc_bike_in['orog'][:,:]
        thk = nc_bike_in['lithk'][:,:]
        frac = nc_bike_in['sftgif'][:,:]
        calv = nc_bike_in['licalvf'][:,:]
        grid = Uniform2DGrid(nc_bike_in['x'][:].data, 
                            nc_bike_in['y'][:].data)
        mask = np.where(frac > 1e-4, True, False)

        gr_mask = None

    elif cfg['icesheet'] == 'AIS' and "cx209" not in nc_file:

        topo = nc_bike_in['orog'][:,:]
        thk = nc_bike_in['lithk'][:,:]
        frac = nc_bike_in['sftgif'][:,:]
        calv = nc_bike_in['licalvf'][:,:]
        grid = Uniform2DGrid(nc_bike_in['x'][:].data, 
                             nc_bike_in['y'][:].data)
        mask = np.where(frac > 1e-4, True, False)

        gr_mask = None

    # for TerraFIRMA ism_nc_file
    elif cfg['icesheet'] == 'AIS' and "cx209" in nc_file:
        
        topo = nc_bike_in['Z_surface'][:,:]
        thk = nc_bike_in['thickness'][:,:]
        frac = nc_bike_in['thickness'][:,:]
        frac = np.where(frac > 1e-4, 1.0, 0.0)
        calv = nc_bike_in['calvingFlux'][:,:]
        grid = Uniform2DGrid(nc_bike_in['x'][:].data, 
                             nc_bike_in['y'][:].data)
        mask = np.where(frac > 1e-4, True, False)

        gr_mask = nc_bike_in['dragCoef'][:,:]

        # In the bisicles*-plot.nc file used here:
        # Ice shelf dragCoef = 0, grounded ice dragCoef > 1
        # non-ice (open ocean) dragCoef = 1

        gr_mask = np.where(gr_mask > 1.00001, True, False)

    else:
        raise ValueError(f"Unsupported ice sheet and nc_file combination: {cfg['icesheet']}, {nc_file}")

    return grid, topo, thk, calv, frac, mask, gr_mask

# Read atmosphere grid and variables
def read_atm_nc(nc_file):
   
    nc_um_in = Dataset(nc_file,'r')

    grid = Uniform2DGrid(nc_um_in['longitude'][:].data,
                            nc_um_in['latitude'][:].data)

    def prep_um(arr):
        return arr[:,:,:]
    
    area = prep_um(nc_um_in['tile_surface_area'])
    spy = 3600*360*24
    smb = prep_um(nc_um_in['ice_smb']) / 918 * spy # values suggest kg/s
    stemp = prep_um(nc_um_in['ice_stemp']) + 273.15
    snow = prep_um(nc_um_in['nonice_snowdepth'])
    shflx = prep_um(nc_um_in['snow_ice_hflux'])
    z_id = nc_um_in['tile_id']

    lon, lat = grid.coords
    nlat, nlon = lat.shape
    nec = z_id.shape[0]
    
    
    # These are passed to Glint from command line in uncicles
    topo_mid = np.array([31.96, 301.57, 550.11, 850.85, 1156.89,
                              1453.16, 1806.74, 2268.25, 2753.70, 3341.4])

    # These are the values that are hardcoded into Glint for nec=10
    topo_max = np.array([0.0, 200.0, 400.0, 700.0, 1000.0, 1300.0, 1600.0, 2000.0, 2500.0, 3000.0, 10000.0])
    
    topo = np.zeros(smb.shape)
    for ec in range(0,nec):
        topo[ec,:,:] = topo_mid[ec]  
        
    return grid, area, topo, topo_max, smb, stemp, snow, shflx

# Function to get files to process
def get_files_to_process(input_dir, start_pattern, num_files=None):

    all_files = sorted(os.listdir(input_dir))
    all_files = [file for file in all_files if file.startswith(start_pattern)]
    all_files_abs = [os.path.join(input_dir, file) for file in all_files]
    if num_files is not None:
        all_files_abs = all_files_abs[:num_files]
    return all_files_abs


# test_TerraFIRMA main function
def test_terrafirma(cfg):

    if cfg['icesheet'] == 'AIS':
        up_tr, down_tr = up_down_pair(PROJ_ANTARCTIC_3031)

    elif cfg['icesheet'] == 'GrIS':
        up_tr, down_tr = up_down_pair(PROJ_ARCTIC_4326)

    else:
        raise ValueError(f"Unsupported ice sheet: {cfg['icesheet']}")

    grid_ism, topo_ism, thk_ism, calv_ism, frac_ism, mask_ism, gr_mask_ism = \
        read_ism_nc(cfg)
    
    input_file_path = f"{cfg['input_dir']}/{cfg['suite_id']}/"

    files = get_files_to_process(input_file_path, start_pattern='atmos', num_files=cfg.get('num_files', None))

    if cfg['gr_fl_mask']:
        total_smb = np.ndarray(shape=(len(files),4))
    
    elif not cfg['gr_fl_mask']:
        total_smb = np.ndarray(shape=(len(files),2))

    else:
        raise ValueError(f"Unsupported gr_mask flag: {cfg['gr_fl_mask']}")
    
    count = 0

    for file in files:
        
        print(f"Processing file: {file}")
        print(f"File {count+1}/{len(files)}")

        grid_um, area_um, topo_um, topo_max, smb_um, stemp_um, snow_um, shflx_um,  = \
            read_atm_nc(file)

        delta_t = 1.0 # typical?

        smb, stemp, delta_snow, shflx \
        = atm_to_ism(smb_um, stemp_um, snow_um, shflx_um,
                     topo_um, area_um, grid_um,
                     topo_ism, thk_ism, frac_ism, mask_ism, grid_ism,
                     up_tr, down_tr, delta_t)
        
        smb += delta_snow / delta_t

        # Convert from m/a to Gt/a in the same way as for TerraFIRMA analysis

        if cfg['icesheet'] == 'GrIS':
            dx = 4800

        elif cfg['icesheet'] == 'AIS':
            dx = 8000

        else:
            raise ValueError(f"Unsupported ice sheet: {cfg['icesheet']}")
        
        smb = smb*(dx**2)*918*1e-12


        total_smb[count,0] = count
        total_smb[count,1] = np.sum(smb)

        if cfg['gr_fl_mask']:

            total_smb[count,2] = ma.sum(ma.masked_array(smb, ~gr_mask_ism))
            total_smb[count,3] = ma.sum(ma.masked_array(smb, gr_mask_ism))
        
        print(f"Total SMB for file {count+1}/{len(files)}: {total_smb[count,1]:.4f} Gt/a")

        count += 1

    return total_smb

def save_results(cfg, smb_array, output_dir):

    time = smb_array[:,0]

    if cfg["gr_fl_mask"]:
        ds = xr.Dataset(
        {
            "total_smb": ("time", smb_array[:,1]),
            "total_smb_grounded": ("time", smb_array[:,2]),
            "total_smb_floating": ("time", smb_array[:,3])
        },
        coords = {"time": time} 
        )
        
    elif not cfg["gr_fl_mask"]:
        ds = xr.Dataset(
            {
                "total_smb": ("time", smb_array[:,1])
            },
            coords = {"time": time}
        )

    else:
        raise ValueError(f"Unsupported gr_mask flag: {cfg['gr_fl_mask']}")

    ds.to_netcdf(os.path.join(output_dir, 'pyglint_smb_data.nc'))

# Main test function
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, help='Path to YAML config')
    args = parser.parse_args()

    # Load config
    cfg = load_config(args.config)

    # Label output dir with branch, config and timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    branch_name = get_git_info()['git_branch']
    config_name = os.path.splitext(os.path.basename(args.config))[0]
    output_dir = os.path.join('test_outputs', f"{branch_name}_{config_name}_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)

    # Redirect stdout/stderr to a log file while still printing to screen
    log_path = os.path.join(output_dir, 'log.txt')
    class Logger(object):
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "w")
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    sys.stdout = sys.stderr = Logger(log_path)

    # Save config copy
    with open(os.path.join(output_dir, 'config_used.yaml'), 'w') as f:
        yaml.dump(cfg, f)

    # Save Git metadata
    git_info = get_git_info()
    with open(os.path.join(output_dir, 'metadata.json'), 'w') as f:
        json.dump(git_info, f, indent=2)

    # Run the test
    print("Running test_TerraFIRMA test with config:")
    print(cfg)
    total_smb = test_terrafirma(cfg)

    # Save results
    save_results(cfg, total_smb, output_dir)

if __name__ == "__main__":
    main()


