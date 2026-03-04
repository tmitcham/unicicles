from netCDF4 import Dataset
from transformation import up_down_pair, Uniform2DGrid, PROJ_ANTARCTIC_3031

bike_nc_file = "bike_cx209c_18510101_plot-AIS-lev0-epsg3031.nc"

nc_bike_in = Dataset(bike_nc_file,'r')

topo = nc_bike_in['Z_surface'][:,:]
x = nc_bike_in['x'][:].data
y = nc_bike_in['y'][:].data

grid = Uniform2DGrid(x, y)

up_tr, down_tr = up_down_pair(PROJ_ANTARCTIC_3031)

lon, lat = up_tr(*grid.coords)

# Write out x/y, lat/lon and topo to new netCDF file

out_nc_file = "bisicles_AIS_cx209_18510101_lat_lon_Z_surface.nc"

with Dataset(out_nc_file, 'w', format='NETCDF4') as nc_out:
    
    # Create dimensions
    nc_out.createDimension('x', len(x))
    nc_out.createDimension('y', len(y))

    # Create variables
    x_var = nc_out.createVariable('x', 'f4', ('x',))
    y_var = nc_out.createVariable('y', 'f4', ('y',))
    lat_var = nc_out.createVariable('lat', 'f4', ('y','x'))
    lon_var = nc_out.createVariable('lon', 'f4', ('y','x'))
    topo_var = nc_out.createVariable('Z_surface', 'f4', ('y','x'))

    # Write data
    x_var[:] = x
    y_var[:] = y
    lat_var[:,:] = lat
    lon_var[:,:] = lon
    topo_var[:,:] = topo