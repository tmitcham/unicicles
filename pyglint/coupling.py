#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 16:10:16 2025

@author: ggslc
"""
import numpy as np
from downscale import interp_to_local, interp_to_surface
from upscale import mean_to_global_mec
from transformation import local_to_global_map, cell_id, fraction_covered
from transformation import  Uniform2DGrid, Box, GlobalLocalGridPair


def delta_h_snow(h_ice, h_snow, thresh=50.0):
    """Calculate change in snow depth due to snow-to-ice conversion.
    
    In thin ice areas (below threshold), convert thick snow to ice
    or convert all ice to snow where snow is thin (below threshold). 

    How it was in GLINT even if it seems a bit odd.

    Parameters
    ----------
    h_ice : ndarray
        Ice thickness (m).
    h_snow : ndarray
        Snow depth (m).
    thresh : float, optional
        Ice thickness threshold (m) below which snow converts to ice. Default is 50.0 m.
    
    Returns
    -------
    ndarray
        Change in snow depth (m). Negative values indicate snow loss due to conversion.
    """
    return np.where( h_ice < thresh,  
                      np.where ( h_snow > thresh, - h_snow, h_ice), 0.0)
    
def glint_conservation_adjust(f2d_ism_nc, mask_ism,  
                              f3d_atm, area_atm,  valid_atm,
                              grid_pair):
    """Enforce conservation of atmospheric fields when downscaling to ice sheet grid.
    
    Adjusts local ice-sheet grid values so that the total field within each 
    atmosphere grid cell equals the vertically-integrated atmospheric column total.
    
    Parameters
    ----------
    f2d_ism_nc : ndarray
        2D field on ice sheet grid (e.g., surface mass balance).
    mask_ism : ndarray
        Boolean mask indicating ice-covered cells on ice sheet grid.
    f3d_atm : ndarray
        3D field on atmospheric grid, shape (lev, N, M).
    area_atm : ndarray
        Cell areas on atmospheric grid.
    valid_atm : ndarray
        Boolean mask of valid atmospheric data, shape (lev, N, M).
    grid_pair : GlobalLocalGridPair
        Mapping between atmospheric and ice sheet grids.
    
    Returns
    -------
    np.ma.MaskedArray
        Adjusted 2D field on ice sheet grid with conservation applied.
    dict
        conservation statistics
    """
    
    # conservation enforcement - see glint_downscaling_gcm
    # roughly - adjust smb so that the total smb_ism within an
    # atmosphere grid box equals the smb_atm integrated vertically.
    # this resembles the glint code. blocky.
    # there is a NASTY! comment in the glint code
    # optmized with the help of co-pilot  

    conservation_stats = {}

    f2d_ism = np.array(f2d_ism_nc, copy=True)

    ism_to_atm_map = grid_pair.local_to_global_map
    grid_ism = grid_pair.local_grid
    grid_atm = grid_pair.global_grid

    dx, dy = grid_ism.spacing 
    area_factor = dx * dy

    # Number of atmosphere cells expected (cell_id produces indices in [0, N*M-1])
    N, M = grid_atm.array_shape
    atm_cell_count = N * M

    # Compute vertical-integrated atm column totals safely: zero out invalid vertical entries
    # f3d_atm and area_atm shapes: (lev, N, M) -> sum over axis=0 results in (N, M)
    atm_col_total_2d = np.sum(np.where(valid_atm, f3d_atm * area_atm, 0.0), axis=0)

    conservation_stats['atm_total'] = np.sum(atm_col_total_2d)

    # Flatten to length N*M with C-order so that index cc = J*M + I matches cell_id(I,J,M,N)
    atm_col_total_flat = atm_col_total_2d.ravel(order='C')

    # Prepare flattened local mappings and local smb
    mapped_ids = ism_to_atm_map.ravel(order='C').astype(np.int64)
    f2d_flat = f2d_ism.ravel(order='C')
    m2d_flat = np.where(mask_ism.ravel(order='C'),1.0,0.0)

    # Sum f2d_ism per atm-cell and count number of local cells per atm-cell
    ism_total_flat = area_factor * np.bincount(mapped_ids, weights=f2d_flat*m2d_flat, minlength=atm_cell_count)
    ism_area_flat = area_factor * np.bincount(mapped_ids, minlength=atm_cell_count)
    ism_ice_area_flat = area_factor * np.bincount(mapped_ids, weights=m2d_flat, minlength=atm_cell_count)

    total_ism = np.sum( f2d_flat * m2d_flat) * area_factor 
    total_atm = np.sum(atm_col_total_flat, where=(ism_ice_area_flat > 0.0))
 
    conservation_stats['ism_total_input'] = total_ism
    conservation_stats['atm_total_masked_to_ism'] = total_atm

    # Compute adjustments only for atm cells that have any valid vertical data and have non-zero ism area
    atm_valid_flat = np.any(valid_atm, axis=0).ravel(order='C')

    delta_flat = np.zeros(atm_cell_count, dtype=float)
    nonz = atm_valid_flat & (ism_ice_area_flat > 0.0)

    if np.any(nonz):
        atm_total_over_ism = np.sum(atm_col_total_flat[nonz])

        delta_flat[nonz] = (atm_col_total_flat[nonz] *  \
                                ism_ice_area_flat[nonz] / ism_area_flat[nonz]    \
                                    - ism_total_flat[nonz]) / ism_ice_area_flat[nonz]

    # Apply adjustments to local grid: each local cell gets the delta of its mapped atm cell
    # For unmapped local cells (flat_map < 0) we do not apply any adjustment (delta_contrib = 0)
    delta_per_local = np.zeros_like(f2d_flat, dtype=float)
    # Only assign for mapped local cells to avoid indexing errors
    delta_per_local = delta_flat[mapped_ids]

    # Update flattened f2d and reshape back
    f2d_flat += delta_per_local
    
    total_ism = np.sum( f2d_flat * m2d_flat) * area_factor 
    conservation_stats['ism_total_after_atm_cell_conservation'] = total_ism
    
    if total_ism > 1e-12:

        f2d_flat*=total_atm/total_ism
    
    total_ism = np.sum( f2d_flat * m2d_flat) * area_factor 
    conservation_stats['ism_total_after_atm_region_conservation'] = total_ism

    f2d_ism = f2d_flat.reshape(f2d_ism.shape, order='C')

    return np.ma.masked_array(f2d_ism, ~mask_ism), conservation_stats



def crop_global(arrs_global, grid_atm, grid_ism, up_transform):
    """Crop global atmospheric fields to region containing ice sheet domain.
    
    Extracts the minimal atmospheric region needed to cover the ice sheet
    domain, reducing computational cost for downscaling operations.
    
    Parameters
    ----------
    arrs_global : list of ndarray
        Global atmospheric fields to crop, all 3D arrays ordered (z,y,x)
    grid_atm : Uniform2DGrid
        Atmospheric grid description.
    grid_ism : Uniform2DGrid
        Ice sheet grid description.
    up_transform : callable
        Transformation function from ice sheet to atmospheric coordinates.
    
    Returns
    -------
    list
        Cropped fields, cropped atmospheric grid, and crop indices (ilo, ihi, jlo, jhi).
    """
    
    # subset
    lon_ism, lat_ism = up_transform(*grid_ism.coords)
    dlon, dlat = [np.max(np.abs(ll[1:] - ll[:-1])) for ll in  grid_atm.axes]
    bb = Box( [np.min(lon_ism) - dlon, np.min(lat_ism) - dlat],
             [np.max(lon_ism) + dlon, np.max(lat_ism) + dlat])

    #latitude crop
    ilo = max(0,np.argmin( np.abs(grid_atm.axes[1] - bb.lo[1])) - 1)
    ihi = min(np.argmin( np.abs(grid_atm.axes[1] - bb.hi[1])) + 1,grid_atm.axes[1].shape[0] )
    
    #longitude crop
    #jlo = max(0,np.argmin( np.abs(grid_atm.axes[0] - bb.lo[0])) - 1)
    #jhi = min(np.argmin( np.abs(grid_atm.axes[0] - bb.hi[0])) + 1,grid_atm.axes[0].shape[0] )
    jlo = 0
    jhi = grid_atm.axes[0].shape[0]

    grid_atm_sub = Uniform2DGrid( grid_atm.axes[0][jlo:jhi], grid_atm.axes[1][ilo:ihi])
    #print (f'cropping {grid_atm.array_shape} -> {grid_atm_sub.array_shape}')

    arrs = [arr[:,ilo:ihi,jlo:jhi] for arr in arrs_global]
    arrs.append(grid_atm_sub)
    arrs.append(ilo)
    arrs.append(ihi)

    return arrs





def atm_to_ism(smb_atm, stemp_atm, snow_atm, shflx_atm,  
               topo_atm, area_atm,  grid_atm,
               topo_ism, lithk_ism, frac_ism, mask_ism, grid_ism,
               up_transform,  down_transform, time_step, 
               conservation=True, snow_to_ice=True, verbose=True): 
    """Downscale atmospheric fields to ice sheet grid.
    
    Interpolates atmospheric surface mass balance, temperature, snow, and heat flux
    to the ice sheet grid.

    Optionally applies conservation adjustment and snow-to-ice conversion.

    Parameters
    ----------
    smb_atm : np.ma.MaskedArray
        Surface mass balance on atmospheric grid (kg m-2 s-1).
    stemp_atm : np.ma.MaskedArray
        Surface temperature on atmospheric grid (K).
    snow_atm : np.ma.MaskedArray
        Snow mass on atmospheric grid (kg m-2).
    shflx_atm : np.ma.MaskedArray
        Sensible heat flux on atmospheric grid (W m-2).
    topo_atm : ndarray
        Topography on atmospheric grid (m).
    area_atm : np.ma.MaskedArray
        Cell areas on atmospheric grid (m2).
    grid_atm : Uniform2DGrid
        Atmospheric grid description.
    topo_ism : ndarray
        Topography on ice sheet grid (m).
    lithk_ism : ndarray
        Ice thickness on ice sheet grid (m).
    frac_ism : np.ma.MaskedArray
        Ice fraction on ice sheet grid.
    mask_ism : ndarray
        Boolean mask of ice-covered cells on ice sheet grid.
    grid_ism : Uniform2DGrid
        Ice sheet grid description.
    up_transform : callable
        Transformation from ice sheet to atmospheric coordinates.
    down_transform : callable
        Transformation from atmospheric to ice sheet coordinates.
    time_step : float
        Time step for integration (s).
    conservation : bool, optional
        Whether to apply conservation adjustment. Default is True.
    snow_to_ice : bool, optional
        Whether to convert snow to ice in thin ice areas. Default is True.

    Returns
    -------
    tuple of np.ma.MaskedArray
        (smb_ism, stemp_ism, delta_snow_ism, shflx_ism) - downscaled fields on ice sheet grid.
    """
    
    smb_atm, stemp_atm, snow_atm, shflx_atm, topo_atm, area_atm, grid_atm, \
    ilo, ihi = \
        crop_global((smb_atm, stemp_atm, snow_atm, shflx_atm, topo_atm, area_atm),
                    grid_atm, grid_ism, up_transform)

    grid_pair = GlobalLocalGridPair(grid_atm, grid_ism, up_transform, down_transform)  

    atm_coords = down_transform(*grid_atm.coords)

    valid_atm = area_atm.data > 0.0
    smb_xyz, snow_xyz, shflx_xyz  = [interp_to_local(f_atm.data, atm_coords, grid_ism.coords,
                              order=1, valid=valid_atm & ~f_atm.mask) 
                         for f_atm in [smb_atm, snow_atm, shflx_atm]]

    stemp_xyz = interp_to_local(stemp_atm.data, atm_coords, grid_ism.coords,
                               order=1, valid=~stemp_atm.mask)

    # this is probably uniform in x,y , but ...
    topo_xyz = interp_to_local(topo_atm, atm_coords, grid_ism.coords, order=1)


    stemp_ism, smb_ism, snow_ism, shflx_ism = \
        [np.ma.masked_array(interp_to_surface(f_xyz, topo_xyz, topo_ism), ~mask_ism)
         for f_xyz in [stemp_xyz, smb_xyz, snow_xyz, shflx_xyz]]

    if conservation and not isinstance(area_atm, type(None)):
        (smb_ism, smb_stats), (snow_ism, _), (shflx_ism, _) = [glint_conservation_adjust(
                            f_ism, mask_ism, 
                            f_atm, area_atm, valid_atm,
                            grid_pair) 
                for f_ism, f_atm in zip([smb_ism, snow_ism, shflx_ism],
                                 [smb_atm, snow_atm, shflx_atm],)]

        if verbose:
            print('SMB conservation stats')
            for key, val in smb_stats.items():
                print (f'{key} =  {smb_stats[key]}')      
        
    if snow_to_ice:  
        #snow to ice conversion - store *change in* snow depth in snow_ism    
        delta_snow_ism = np.ma.masked_array(delta_h_snow(lithk_ism, snow_ism), ~mask_ism)
        smb_ism -= delta_snow_ism / time_step 
    else:
        delta_snow_ism = np.ma.masked_array( np.zeros_like(smb_ism.data), ~mask_ism)

    return smb_ism, stemp_ism, delta_snow_ism, shflx_ism


def splice_global(global_arr, region_arr, frac_cover,  ilo, ihi):
    """Blend regional ice sheet data into global atmospheric field.
    
    Replaces global field values in region with weighted blend of regional
    and global data, using ice coverage fraction as weight.
    
    Parameters
    ----------
    global_arr : ndarray
        Global field array, 3D ordered(z,lat,lon)
    region_arr : ndarray
        Regional field array, 3D ordered (z,lat,lon).
    frac_cover : ndarray
        Ice coverage fraction, 2D ordered(lat, lon) 
    ilo : int
        Lower latitude index of region in global array.
    ihi : int
        Upper latitude index of region in global array.
    
    Returns
    -------
    ndarray
        Modified global array with regional data blended in.
    """
    
    global_arr[:,ilo:ihi,:] = region_arr*frac_cover + (1.0-frac_cover)*global_arr[:,ilo:ihi,:]

    return global_arr

def ism_to_atm(topo_max_atm, area_atm, grid_atm,
               ice_frac_atm, delta_snow_atm, hflux_atm, calv_atm, topo_atm, 
               topo_ism, frac_ism, 
               hflx_ism, d_snow_ism, calv_ism,
               mask_ism, grid_ism,
               up_transform, down_transform):
    """Upscale ice sheet fields to atmospheric grid.
    
    Computes area-weighted means of ice sheet fields (ice fraction, topography,
    snow depth, heat flux, calving) on the atmospheric grid.
    
    Parameters
    ----------
    topo_max_atm : ndarray
        Maximum topography on atmospheric grid (m).
    area_atm : np.ma.MaskedArray
        Cell areas on atmospheric grid (m2).
    grid_atm : Uniform2DGrid
        Atmospheric grid description.
    ice_frac_atm : np.ma.MaskedArray
        Ice fraction on atmospheric grid (output array).
    delta_snow_atm : np.ma.MaskedArray
        Change in snow depth on atmospheric grid (m, output array).
    hflux_atm : np.ma.MaskedArray
        Heat flux on atmospheric grid (W m-2, output array).
    calv_atm : np.ma.MaskedArray
        Calving flux on atmospheric grid (kg m-2 s-1, output array).
    topo_atm : np.ma.MaskedArray
        Topography on atmospheric grid (m, output array).
    topo_ism : ndarray
        Ice sheet surface topography (m).
    frac_ism : np.ma.MaskedArray
        Ice coverage fraction on ice sheet grid.
    hflx_ism : np.ma.MaskedArray
        Heat flux on ice sheet grid (W m-2).
    d_snow_ism : np.ma.MaskedArray
        Snow depth change on ice sheet grid (m).
    calv_ism : np.ma.MaskedArray
        Calving flux on ice sheet grid (kg m-2 s-1).
    mask_ism : ndarray
        Boolean mask of ice-covered cells on ice sheet grid.
    grid_ism : Uniform2DGrid
        Ice sheet grid description.
    up_transform : callable
        Transformation from ice sheet to atmospheric coordinates.
    down_transform : callable
        Transformation from atmospheric to ice sheet coordinates.
    
    Returns
    -------
    list
        Updated atmospheric fields: [ice_frac_atm, topo_atm, delta_snow_atm, hflux_atm, calv_atm].
    """

    # GLINT 3D ouputs : gtopo - ice sheet surface,
    #                   gfrac - ice covered fraction
    #                   grofi - ice run off (calving?)
    #                   grofl - liquid run off
    #                   ghflx - heat flux
    #                   gsdep - snow depth (in total / out anomaly)
    #                   gcalv = calving flux
    # GLINT 0D outputs : ice volume

    # wrappers/ukesm-ice_NETCDF/gl_mod.f90 calls glint to obtain
    # gsdep (in also), gfrac, ghflx, gcalv, glfrac, ice_volume
    # i.e ignore grofi , grofl

    dbg = 0

    area_atm_sub, grid_atm_sub, ilo, ihi = \
        crop_global([area_atm], grid_atm, grid_ism, up_transform)
        
    grid_pair = GlobalLocalGridPair(grid_atm_sub, grid_ism, up_transform, down_transform)    
        

    region_data = mean_to_global_mec([frac_ism.data, 
                                      topo_ism.data, 
                                      d_snow_ism.data, 
                                      hflx_ism.data, 
                                      calv_ism.data], 
                                     [True, False, False, False, False], 
                                     topo_ism.data, mask_ism,
                                     topo_max_atm, 
                                     grid_pair.local_to_global_map,
                                     area_atm_sub.shape)
   
    
    frac_cover = fraction_covered(grid_pair)
    global_data = [ice_frac_atm, topo_atm, delta_snow_atm, hflux_atm, calv_atm]
    global_data = [splice_global(g, r , frac_cover, ilo, ihi)
                   for g,r in zip(global_data, region_data)]

    return global_data





