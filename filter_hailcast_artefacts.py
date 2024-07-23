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
path to the input directory
path to the output directory
start day in format YYYYMMDD
end day in format YYYYMMDD

optional:
threshold for hailcast data in mm/h (default 2 mm/h)
kernel size for binary dilatation in km (default 17 km)
number of time steps to consider before and after current time step in time steps (default 4)


---------------------------------------------------------
OUT

---------------------------------------------------------
EXAMPLE CALL
python /home/kbrennan/filter_hailcast_artefacts/filter_hailcast_artefacts.py /net/litho/atmosdyn2/kbrennan/data/climate/present/5min_2D/ /net/litho/atmosdyn2/kbrennan/data/climate/present/filtered_hailcast/ 20210101 20220101
---------------------------------------------------------
Killian P. Brennan
14.05.2024
---------------------------------------------------------
todos:
change data source to the 5min data
add more temporal ambiguity (+/- 4 time steps)
add threshold on slope angle (inflated)

"""

import os
import json
import xarray as xr
import numpy as np
import pandas as pd
from argparse import ArgumentParser

import scipy.ndimage as ndimage
from skimage.morphology import disk


def filter_hailcast_artefacts(
    inpath,
    outpath,
    start_day,
    end_day,
    temporal_ambiguity=4,  # time steps
    threshold=2,  # mm/h
    kernel_size_km=17,  # km
):
    start_day = pd.to_datetime(start_day, format="%Y%m%d")
    end_day = pd.to_datetime(end_day, format="%Y%m%d")

    daylist = pd.date_range(start_day, end_day)
    day_m_1 = start_day - pd.DateOffset(days=1)

    # make output directory
    os.makedirs(outpath, exist_ok=True)

    for i, day in enumerate(daylist):
        today_str = day.strftime("%Y%m%d")
        print(today_str)
        if i == 0:
            yesterday = xr.open_dataset(
                os.path.join(inpath, "lffd" + day_m_1.strftime("%Y%m%d") + "_0606.nz")
            )
        else:
            yesterday = today
        if i == 0:
            today = xr.open_dataset(
                os.path.join(inpath, "lffd" + today_str + "_0606.nz")
            )
        else:
            today = tomorrow

        tomorrow_path = os.path.join(
            inpath, "lffd" + daylist[i + 1].strftime("%Y%m%d") + "_0606.nz"
        )
        if os.path.exists(tomorrow_path):
            tomorrow = xr.open_dataset(tomorrow_path)
        else:
            print('warning: last day in time range, filling with high precip value - no hail will be filtered for the last n=spacial_ambiguity time steps')
            # for the last day, the next day does not exist, fill with high precip value so no hail is filtered
            tomorrow = xr.full_like(today, 100)

        combined = xr.concat(
            [
                yesterday.isel(time=slice(-(temporal_ambiguity + 1), -1)),
                today,
                tomorrow.isel(time=slice(0, temporal_ambiguity + 1)),
            ],
            dim="time",
        )
        TOT_PREC = combined.TOT_PREC
        DHAIL_MX = today.DHAIL_MX
        filtered = filter_daily_field(
            DHAIL_MX,
            TOT_PREC,
            threshold=threshold,
            kernel_size_km=kernel_size_km,
            temporal_ambiguity=temporal_ambiguity,
        )

        # write compressed netcdf
        filtered.to_netcdf(
            os.path.join(outpath, "lffd" + today_str + "_0606.nc"),
            format="NETCDF4",
            engine="netcdf4",
            encoding={"DHAIL_MX": {"zlib": True, "complevel": 9}},
        )
    return


# def filter_5min_climate(inpath, outpath, start_day, end_day):
#     temporal_ambiguity = 4  # time steps
#     threshold = 2  # mm/h
#     kernel_size = 8  # grid points

#     # list all files in the input directory
#     files = os.listdir(inpath)
#     files = [f for f in files if f.endswith(".nc")]
#     files = [f for f in files if f.startswith("lffd")]
#     files = sorted(files)
#     # format is lffdYYYYMMDDHHmmss.nc

#     for i, f in enumerate(files):
#         if int(f[4:12]) >= (int(end_day)):
#             print("reached end of specified time range")
#             break
#         print(f)
#         if int(f[4:12]) < int(start_day):
#             print("skipping, not in specified time range")
#             continue
#         # print("processing")
#         ds = xr.open_mfdataset(
#             [
#                 os.path.join(inpath, f)
#                 for f in files[i - temporal_ambiguity : i + temporal_ambiguity + 1]
#             ]
#         )
#         rr = ds.TOT_PREC.max(dim="time") * 12
#         DHAIL_MX = ds.DHAIL_MX.isel(time=temporal_ambiguity)
#         m = build_mask(rr, threshold=threshold, kernel_size=kernel_size)
#         DHAIL_MX_filtered = DHAIL_MX.where(m, 0)
#         # add back the time dimension
#         # DHAIL_MX_filtered = DHAIL_MX_filtered.expand_dims(dim={"time": ds.DHAIL_MX.isel(time=temporal_ambiguity).time})
#         DHAIL_MX_filtered = DHAIL_MX_filtered.expand_dims("time")
#         DHAIL_MX_filtered.to_netcdf(os.path.join(outpath, f))

#     return


def build_mask(rr, threshold=5, kernel_size=8):
    kernel = disk(kernel_size)
    m = rr > threshold
    m = ndimage.binary_dilation(m, structure=kernel)
    return m


def filter_daily_field(
    DHAIL_MX,
    TOT_PREC,
    threshold=2,  # mm/h
    kernel_size_km=17,  # km
    temporal_ambiguity=4,  # time steps
):
    """
    Filter out hailcast artefacts using a threshold and a binary dilatation
    Parameters
    ----------
    ds : xarray dataset
    hailcast data
    threshold : float
    threshold for hailcast data
    kernel_size_km : int
    kernel size for binary dilatation
    temporal_ambiguity : int
    number of time steps to consider before and after current time step

    Returns
    -------
    ds : xarray dataset
    hailcast data with artefacts removed

    """
    DHAIL_MX_filtered = xr.full_like(DHAIL_MX, 0)

    # get dt from df.time
    dt = (
        pd.to_timedelta(DHAIL_MX.time.diff(dim="time").mean().values).total_seconds()
        / 3600
    )
    # dt = 1 / 12  # 5 min

    # get grid spacing from df.lat
    r_earth = 6371  # km
    lat = DHAIL_MX.lat.values
    dlat = np.abs(lat[1] - lat[0])
    dx = r_earth * np.radians(dlat)
    # get grid spacing in km
    dx = dx.mean()
    kernel_size_gp = int(np.round(kernel_size_km / dx))

    mask = xr.full_like(DHAIL_MX, False)

    for i in range(len(DHAIL_MX.time)):
        # tot_prec start index is offset from dhail start index by temporal ambiguity
        rr = (
            TOT_PREC.isel(time=slice(i, i + 2 * temporal_ambiguity + 1)).max(dim="time")
            / dt
        )
        m = build_mask(rr, threshold=threshold, kernel_size=kernel_size_gp)
        mask[i] = m.copy()

    DHAIL_MX_filtered = DHAIL_MX.where(mask, 0)

    # add filtering attributes
    DHAIL_MX_filtered.attrs["filter_threshold_mm_per_h"] = threshold
    DHAIL_MX_filtered.attrs["filter_kernel_size_gp"] = kernel_size_gp
    DHAIL_MX_filtered.attrs["filter_kernel_size_km"] = kernel_size_km
    DHAIL_MX_filtered.attrs["filter_temporal_ambiguity_time_steps"] = temporal_ambiguity

    return DHAIL_MX_filtered


if __name__ == "__main__":

    parser = ArgumentParser(description="filter orographic hailcast artefacts")

    parser.add_argument("inpath", type=str, help="path to the input directory")
    parser.add_argument("outpath", type=str, help="path to the output directory")
    parser.add_argument("start_day", type=str, help="start day in format YYYYMMDD")
    parser.add_argument("end_day", type=str, help="end day in format YYYYMMDD")
    parser.add_argument(
        "--threshold", type=float, default=2, help="threshold for hailcast data"
    )
    parser.add_argument(
        "--kernel_size_km",
        type=int,
        default=17,
        help="kernel size for binary dilatation",
    )
    parser.add_argument(
        "--temporal_ambiguity",
        type=int,
        default=4,
        help="number of time steps to consider before and after current time step",
    )

    args = parser.parse_args()

    filter_hailcast_artefacts(
        args.inpath,
        args.outpath,
        args.start_day,
        args.end_day,
        threshold=args.threshold,
        kernel_size_km=args.kernel_size_km,
        temporal_ambiguity=args.temporal_ambiguity,
    )
