#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
---------------------------------------------------------
filter orographic hailcast artefacts

should work for operational mode, as well as for the climate simulations
algorithm working principle
    max over previous and following time steps
    apply threshold
    apply binary dilatation
    apply resulting mask
---------------------------------------------------------
IN

---------------------------------------------------------
OUT

---------------------------------------------------------
EXAMPLE CALL
python /home/kbrennan/filter_hailcast_artefacts/filter_hailcast_artefacts.py /net/litho/atmosdyn2/kbrennan/data/climate/present/5min_2D/ /net/litho/atmosdyn2/kbrennan/data/climate/present/filtered_hailcast/ 20210101 20220101
---------------------------------------------------------
Killian P. Brennan
14.05.2024
---------------------------------------------------------

"""

import os
import json
import xarray as xr
import numpy as np
import pandas as pd
from argparse import ArgumentParser

import scipy.ndimage as ndimage
from skimage.morphology import disk


def filter_hailcast_artefacts(inpath, outpath, start_day, end_day):
    start_day = pd.to_datetime(start_day, format="%Y%m%d")
    end_day = pd.to_datetime(end_day, format="%Y%m%d")

    daylist = pd.date_range(start_day, end_day)

    # make output directory
    os.makedirs(outpath, exist_ok=True)

    for day in daylist:
        day_str = day.strftime("%Y%m%d")
        print(day_str)
        original = xr.open_dataset(os.path.join(inpath,'lffd'+day_str+'_0606.nz'))
        filtered = apply_filter(original)

        filtered.DHAIL_MX.to_netcdf(os.path.join(outpath, 'lffd'+day_str+'_0606.nc'))

    return


def apply_filter(
    ds,
    threshold=5,  # mm/h
    kernel_size=8,  # km
    temporal_ambiguity=2,  # time steps
):
    '''
    Filter out hailcast artefacts using a threshold and a binary dilatation
    Parameters
    ----------
    ds : xarray dataset
    hailcast data
    threshold : float
    threshold for hailcast data
    kernel_size : int
    kernel size for binary dilatation
    temporal_ambiguity : int
    number of time steps to consider before and after current time step

    Returns
    -------
    ds : xarray dataset
    hailcast data with artefacts removed
        
        '''

    ds["DHAIL_MX_filtered"] = ds["DHAIL_MX"]

    # get dt from df.time
    # dt = pd.to_timedelta(ds.time.diff(dim="time").mean().values).total_seconds() / 3600
    dt = 1 / 12  # 5 min

    # get grid spacing from df.lat
    r_earth = 6371  # km
    lat = ds.lat.values
    dlat = np.abs(lat[1] - lat[0])
    dx = r_earth * np.radians(dlat)
    # get grid spacing in km
    dx = dx.mean()
    kernel_size = int(np.round(kernel_size / dx))
    kernel = disk(kernel_size)

    mask = xr.full_like(ds.TOT_PREC, False)

    for i in range(len(ds.time)):
        start = i - temporal_ambiguity
        end = i + temporal_ambiguity + 1
        # make sure the start and end indices are within the time dimension
        if start < 0:
            start = 0
        if end > len(ds.time):
            end = len(ds.time)
        rr = ds.TOT_PREC.isel(time=slice(start, end)).max(dim="time")
        rr = rr / dt
        m = rr > threshold
        m = ndimage.binary_dilation(m, structure=kernel)
        mask[i] = m.copy()

    ds["DHAIL_MX_filtered"] = ds["DHAIL_MX"].where(mask, 0)

    return ds


if __name__ == "__main__":

    parser = ArgumentParser(description="filter orographic hailcast artefacts")

    parser.add_argument(
        'inpath',
        type=str,
        help='path to the input directory'
    )
    parser.add_argument(
        'outpath',
        type=str,
        help='path to the output directory'
    )
    parser.add_argument(
        'start_day',
        type=str,
        help='start day in format YYYYMMDD'
    )
    parser.add_argument(
        'end_day',
        type=str,
        help='end day in format YYYYMMDD'
    )

    args = parser.parse_args()

    filter_hailcast_artefacts(args.inpath, args.outpath, args.start_day, args.end_day)
