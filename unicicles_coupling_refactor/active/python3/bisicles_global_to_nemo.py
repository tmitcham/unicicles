#!/usr/bin/env python

"""
bisicles_global_to_nemo

Carries out two functions, put together here because they both take
regional domain input and produce global fields for NEMO.

--do_iceberg_field
Unifies the regional iceberg calving seed files produced by each instance of
regional_ice_sheet_calving_to_ocean into a global file

--do_shelf_geometry
Reads BISICLES AIS plot file directly and creates global ocean bathymetry
and ice-draft boundary files

Call by suite as task ice_sheets_to_ocean
"""

import sys
import numpy as np
from netCDF4 import Dataset
import cf
from amrfile import io as amrio
import mapping_class as mp
from common_arg_to_file_exist import arg_to_file_exist
from regional_iceberg_seeds_to_global import splice_iceberg_seeds
from common_params_and_constants import bgtn_bathy_threshold, \
    bgtn_max_jindex, bgtn_mindep_bath, bgtn_mincav, bgtn_mindep_isf, \
    bgtn_nflag, bgtn_pflag, bgtn_thresh, bgtn_topopen

THRESH = bgtn_thresh
MINDEP_ISF = bgtn_mindep_isf
MINDEP_BATH = bgtn_mindep_bath
MINCAV = bgtn_mincav


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--do_shelf_geometry",
                        help="calculate new isf_draft for NEMO",
                        action="store_true")
    parser.add_argument("--hdf5_input",
                        help="input, name of bisicles input file")
    parser.add_argument("--nc_output", help="output, name of nemo output file")
    parser.add_argument("--nc_template",
                        help="input, name of nemo original file")
    parser.add_argument("--new_bathy",
                        help="write a new NEMO bathymetry as well as the "
                        "icedraft")
    parser.add_argument("--verbose")
    parser.add_argument("--do_iceberg_field",
                        help="make a global iceberg field for NEMO",
                        action="store_true")
    parser.add_argument("--input_berg_GrIS",
                        help="input, name of regional Greenland ancil")
    parser.add_argument("--input_berg_AIS",
                        help="input, name of regional Antarctic ancil")
    parser.add_argument("--input_berg_global",
                        help="input, name of background global ancil")
    parser.add_argument("--output_berg_global",
                        help="output, name of modified global ancil")
    parser.add_argument("--halo", help="NEMO data contains a halo",
                        action="store_true")
    args = parser.parse_args()

    err = 0

    l_verb = False
    if args.verbose:
        l_verb = True

    l_bath = False
    if args.new_bathy:
        l_bath = True

    do_shelf_geometry = False
    if args.do_shelf_geometry:
        do_shelf_geometry = True

    do_iceberg_field = False
    if args.do_iceberg_field:
        do_iceberg_field = True

    l_halo = False
    if args.halo:
        l_halo = True

    hdf5_input = None
    nc_template = None
    nc_output = None
    if do_shelf_geometry:
        err = 0
        hdf5_input, err = arg_to_file_exist(args.hdf5_input, mandatory=True,
                                            err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)
        nc_template, err = arg_to_file_exist(args.nc_template, mandatory=True,
                                             err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)
        nc_output, err = arg_to_file_exist(args.nc_output, mandatory=True,
                                           io="out", err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)

    GrIS_berg_filename = None
    AIS_berg_filename = None
    global_berg_filename_in = None
    global_berg_filename_out = None
    if do_iceberg_field:
        err = 0
        GrIS_berg_filename, err = arg_to_file_exist(args.input_berg_GrIS,
                                                    mandatory=False, err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)
        AIS_berg_filename, err = arg_to_file_exist(args.input_berg_AIS,
                                                   mandatory=False, err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)
        global_berg_filename_in, err \
            = arg_to_file_exist(args.input_berg_global, err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)
        global_berg_filename_out, err \
            = arg_to_file_exist(args.output_berg_global, io="out", err=err)
        if err != 0:
            parser.print_help()
            sys.exit(err)

    return do_shelf_geometry, hdf5_input, nc_template, nc_output, \
        l_verb, l_bath, do_iceberg_field, GrIS_berg_filename, \
        AIS_berg_filename, global_berg_filename_in, \
        global_berg_filename_out, l_halo


def enforce_minima(array, thresh, mindep):

    #where depth is less than the threshold for existence, set to 0
    zero = np.where((array < thresh))
    array[zero] = 0.

    #where depth is enough to exist but below the practical minimum,
    #set to the practical minimum
    shallow = np.where((array < mindep) & (array > thresh))
    array[shallow] = mindep

    return array


def make_surface_mask(bathy, isf, mindep_isf, mindep_bath, thresh):

    mask = np.zeros_like(bathy)

    #flag the surface as non-ocean where significant ice shelf is present
    #or where the bathymetry is not deep enough for either ocean or ice shelf
    # mask[ (isf > mindep_isf-thresh) | (bathy < mindep_bath+thresh) ] = 1
    mask[(isf > thresh) | (bathy < thresh)] = 1

    return mask


def isf_from_bisicles_median(bathy_in, isf_in, n_mapping, x_mapping,
                             y_mapping, bathy_thresh, l_verb=False):

    """look at the collection of BISICLES points that map to a NEMO cell and
       derive an average bathymetry and water column depth from them, including
       whether there is ice there and if it is grounded """

    n_contrib = 0.
    n_cav = 0
    n_open = 0
    bathy_acc = 0.
    isf_acc = 0.

    #loop over each BISICLES point that maps to this NEMO cell
    for n in range(int(n_mapping)):

        #get the ocean and ice depths from BISICLES
        bathy_in_point = bathy_in[int(y_mapping[n]), int(x_mapping[n])]
        isf_in_point = isf_in[int(y_mapping[n]), int(x_mapping[n])]

        if l_verb:
            print(n, bathy_in_point, isf_in_point)

        if bathy_in_point > bathy_thresh:
            #if the bathymetry is below sea-level, add information from this BISICLES point
            #to the contributions to NEMO
            n_contrib = n_contrib + 1.
            bathy_acc = bathy_acc + bathy_in_point
            isf_acc = isf_acc + isf_in_point

            if (np.abs(isf_in_point - bathy_in_point) > THRESH) \
               & (np.abs(isf_in_point) > THRESH):
                #if BISICLES has significant ice depth and a water column underneath
                #the ice, consider this a point with significant ocean cavity
                #Not grounded, not open ocean
                n_cav = n_cav + 1

            if np.abs(isf_in_point) < THRESH:
                #if BISICLES does not have significant ice depth, consider this as
                # Open ocean
                n_open = n_open + 1

    if l_verb:
        print(n_contrib, n_cav, n_open)
    if l_verb:
        print(isf_acc, bathy_acc)

    isf_out = bgtn_pflag  # Use for detection of out-of-domain points later
    bathy_out = bgtn_pflag

    #if enough BISICLES points in this mapping to give some information to NEMO, work
    #out what information we're going to give
    if n_contrib > 0:

        bathy_out = bathy_acc / n_contrib

        if (n_open + n_cav) >= n_contrib / 2.:
            # Enough points have some kind of ocean that we need to have an ocean column
            # if n_open > n_cav:
            if n_open == n_contrib:
                isf_out = 0.
            else:
                # Enough points are ungrounded and not open ocean, we have a cavity
                # Set water column to the average column depth of all the BISICLES
                # points that map to this NEMO cell
                isf_out = np.abs((bathy_acc - isf_acc)) / n_contrib
        else:
            # Not enough ocean points. Ground ice
            isf_out = bgtn_nflag

        if l_verb:
            print(isf_acc / n_contrib, isf_out)

    return isf_out, bathy_out


def splice_field_distance(full_bath, ice_bath, ice_isf, merge_mask=None):

    """ Merge the regional bathymetry field into a global reference, weighting on
        distance from the ice shelf (or another specified mask)  """

    from scipy import ndimage

    if merge_mask is None:
        # if we haven't been given a mask to base the weights on, make one
        # that reflects where the ice shelf is
        merge_mask = np.zeros_like(ice_isf)
        is_ice = np.where((ice_isf > THRESH) & (ice_isf < 1e4))
        merge_mask[is_ice] = 1.

    # build weights based on distance from mask. "Distance" is simply done
    # by counting grid cells assuming a uniform grid, this is not a great approximation
    # to the real distance for an ORCA grid.
    # /slower/?
    ice_dist_int = ndimage.distance_transform_edt(1 - merge_mask)
    ice_dist = ice_dist_int.astype(np.float)
    is_ice = np.where(ice_dist < 1.)
    ice_dist[is_ice] = 1.
    # Should we just use the distance as a linear factor? Square it? 
    # Multiply by a factor? How should it drop off
    # ice_weights = 1. / ice_dist  # drops off too fast?
    linear_distance_threshold = 10.  # linear dropoff over specified distance
    ice_weights = np.maximum(
        (1 - (ice_dist - 1.) / linear_distance_threshold), 0)

    # flag where we don't believe the BISICLES bathmetry
    bad_ice_bath = np.where(ice_bath > 1e4)
    ice_weights[bad_ice_bath] = 0.

    ny = np.shape(ice_bath)[0]
    nx = np.shape(ice_bath)[1]

    full_bath_new = np.copy(full_bath)

    #merge the regional bathymetry in to the global field using the weights
    full_bath_new[:ny, :nx] = full_bath[:ny, :nx] * \
        (1 - ice_weights) + ice_bath * ice_weights * -1

    full_bath_new = np.maximum(full_bath_new, 0.)

    return full_bath_new


if __name__ == "__main__":

    # Read the command line arguments
    do_shelf_geometry, hdfplotfile, nemo_templatefile, rebin_ncfile, \
        l_verb, l_bath, do_iceberg_field, GrIS_berg_filename, \
        AIS_berg_filename, global_berg_filename_in, \
        global_berg_filename_out, l_halo = parse_commandline()

    if do_iceberg_field:
        # combine the GrIS and AIS iceberg fields (where available) into a global one 
        final_berg = splice_iceberg_seeds(GrIS_berg_filename,
                                          AIS_berg_filename,
                                          global_berg_filename_in)

        cf.write(final_berg, global_berg_filename_out)

    if do_shelf_geometry:
        # nemo_templatefile="BIKE_geom_rebin_eORCA025_grid.nc_match_12.0_30.nc"
        bisicles_to_nemo_mapping_file = "bikegridmap_file.map2d"

        bathy_threshold = bgtn_bathy_threshold
        isf_threshold = bgtn_nflag

        #these are parameters for the uniform grid of information we'll get from BISICLES
        #in particular, "level" needs to match what was used in producing the bisicles_to_nemo_mapping_file
        #that's being used here
        level = 2  # Level of grid refinement
        order = 1  # Interpolation order, 0 for piecewise constant, 1 for linear

        if l_verb:
            print('Reading hdf file')

        #load the ice depth and bathymetry fields from the BISICLES plot file
        amrID = amrio.load(hdfplotfile)
        lo, hi = amrio.queryDomainCorners(amrID, level)
        _, _, isf = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
        _, _, bathy = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)

        # Somewhere in the tangle of reading the hdf5 this way gets me arrays
        # transposed compared to the way they were read/used when generating
        # the original lat and lon metadata.
        bathy = np.transpose(bathy)
        isf = np.transpose(isf)
        ##

        if l_verb:
            print('Restoring mapping files')

        #load the precomputed mapping for how BISICLES grid centres fit into
        #NEMO grid boundaries
        bn_map = mp.load(bisicles_to_nemo_mapping_file)

        #grab the array dimensions we're working with
        ny_n = np.shape(bn_map.x)[0]
        nx_n = np.shape(bn_map.x)[1]

        isf_cav = np.zeros([ny_n, nx_n])
        bathy_ism = np.zeros([ny_n, nx_n])

        if l_verb:
            print("remap geometry to NEMO")

        #loop over each NEMO grid cell
        for j in range(ny_n):
            if l_verb:
                print('\r', j, '  ', end=' ')
            sys.stdout.flush()
            for i in range(nx_n):

                #find the ice iand ocean depths at the BISICLES points that map into 
                #this NEMO cell and produce an ice depth and bathymetry for it
                isf_cav[j, i], bathy_ism[j, i] \
                    = isf_from_bisicles_median(bathy, isf, bn_map.nmap[j, i],
                                               bn_map.x[j, i, :],
                                               bn_map.y[j, i, :],
                                               bathy_threshold, False)

                # No/bad data in mapping
                if isf_cav[j, i] > bgtn_pflag - 1:
                    # Default open ocean
                    isf_cav[j, i] = 0.
                    # 300 is representative of ORCA025, but will work for
                    # ORCA1 too
                    if j < bgtn_max_jindex:
                        # This far south, any bad data is eORCA grid weirdness
                        # shouldn't matter, but nicer if we ground it to
                        # whatever the bathymetry array has
                        isf_cav[j, i] = bgtn_nflag

        # Load the reference bathymetry and ice depth for NEMO
        template = Dataset(nemo_templatefile, 'r')
        bathy_tem_v = template.variables['Bathymetry_isf']
        isf_tem_v = template.variables['isf_draft']
        bathy_tem = bathy_tem_v[:]
        isf_tem = isf_tem_v[:]

        #grab the array dimensions we're working with now
        ny_f = np.shape(isf_tem)[0]
        nx_f = np.shape(isf_tem)[1]

        if l_bath:
            # merge the remapped BISICLES bathymetry info into the global
            # field from the NEMO reference
            bathy = splice_field_distance(l_halo, bathy_tem, bathy_ism, isf_cav)
            # Wrap for halo
            if l_halo:
                bathy[:, -2:] = bathy[:, 0:2]
        else:
            #just use the NEMO reference bathymetry
            bathy = bathy_tem
            # Wrap for halo
            if l_halo:
                bathy[:, -2:] = bathy[:, 0:2]

        # Default to an ice draft field with open ocean, no ice depth
        isf_new = np.zeros([ny_f, nx_f])

        # Do non-open, non-grounded points. Set ice draft to be the bathymetry
        # minus the water column depth that we remapped from BISICLES
        update = np.where(isf_cav > THRESH)
        isf_new[0:ny_n, 0:nx_n][update] = bathy[0:ny_n, 0:nx_n][update] - \
                                          isf_cav[update]

        # Ground ice where the flag says there is no ocean column
        grounded = np.where(isf_cav < 0)
        isf_new[0:ny_n, 0:nx_n][grounded] = bathy[0:ny_n, 0:nx_n][grounded]

        # Wrap for halo
        if l_halo:
            isf_new[:, -2:] = isf_new[:, 0:2]
        # Sanity, safety limit
        isf_new = np.maximum(isf_new, 0)

        isf = isf_new

        # Safety limits for pre-set minimum valid depths
        if l_bath:
            bathy = enforce_minima(bathy, THRESH, MINDEP_BATH)

        isf = enforce_minima(isf, THRESH, MINDEP_ISF)

        # Need to enforce the reference surface ocean/not-ocean mask on the new
        # geometry we've just created. Derive what these are
        mask_tem = make_surface_mask(bathy_tem, isf_tem, MINDEP_ISF,
                                     MINDEP_BATH, THRESH)
        mask = make_surface_mask(bathy, isf, MINDEP_ISF, MINDEP_BATH, THRESH)

        # where do the masks disagree?
        fillmask = mask - mask_tem

        # Reference data is open, new isn't; remove shelf, possibly deepen bathy
        closed = np.where((fillmask > THRESH))
        if l_verb:
            print(np.shape(closed)[1], "were open ocean, now closed - opening")
        isf[closed] = 0

        closed = np.where((fillmask > THRESH) & (bathy < MINDEP_BATH))
        if l_verb:
            print(np.shape(closed)[1], "of which now have land, not ice")
        # if l_bath: bathy[closed] = MINDEP_BATH
        if l_bath:
            bathy[closed] = bgtn_topopen

        # New data is open, reference isn't, bathymetry is shallow ; create shelf
        # and ground it
        open_shallow = np.where( (fillmask < -THRESH) &
                                 (bathy < MINDEP_ISF + MINCAV) )
        if l_verb:
            print(np.shape(open_shallow)[1], "were closed ocean, now open ",
                  "and shallow - creating grounded ice")
        isf[open_shallow] = bathy[open_shallow]

        # New data is open, reference isn't, bathymetry is deep ; create shelf
        open_deep = np.where( (fillmask < -THRESH) &
                              (bathy >= MINDEP_ISF + MINCAV) )
        if l_verb:
            print(np.shape(open_deep)[1], "were closed ocean, now open ",
                  "and deep - creating a cavity")
        # isf[open_deep] = MINDEP_ISF
        isf[open_deep] = bgtn_topopen

        # New data has too thin cavity; ground it
        thin_cavity = np.where( (isf > THRESH) &
                                (np.abs(bathy - isf) > THRESH) &
                                (np.abs(bathy - isf) < MINCAV) )
        if l_verb:
            print(np.shape(thin_cavity)[1],
                  "have a cavity, but I'm grounding as it's thinner than ",
                  MINCAV, " m")
        isf[thin_cavity] = bathy[thin_cavity]

        # Wrap the halo
        if l_halo:
            bathy[:, -2:] = bathy[:, 0:2]
            isf[:, -2:] = isf[:, 0:2]

        # Write out to netCDF file
        output = Dataset(rebin_ncfile, 'w', format='NETCDF3_CLASSIC')
        for dimname in template.dimensions:
            dim = template.dimensions[dimname]

            output.createDimension(dimname, len(dim))

        bathyout = output.createVariable('Bathymetry_isf', bathy_tem_v.dtype,
                                         bathy_tem_v.dimensions)
        for attr in bathy_tem_v.ncattrs():
            bathyout.setncattr(attr, getattr(bathy_tem_v, attr))
        isfout = output.createVariable('isf_draft', isf_tem_v.dtype,
                                       isf_tem_v.dimensions)
        for attr in isf_tem_v.ncattrs():
            isfout.setncattr(attr, getattr(isf_tem_v, attr))
 
        bathyout[:] = bathy[:]
        isfout[:] = isf[:]

        output.close()
