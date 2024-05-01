#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
---------------------------------------------------------
filter orographic hailcast artefacts

should work for operational mode, as well as for the climate simulations
use column integrated graupel (available hourly for climate and weather)
asuming 100km/h maximaum storm propagation velocity, filter out all hailcast signal that is within 40km of graupel signal above threshold
take max of cigraupel from hour before and hour after current timestep
algorithm working principle
    apply threshold
    apply binary dilatation
    apply resulting mask
---------------------------------------------------------
IN

hailcast: xarray dataset
    hailcast data
graupel: xarray dataset
    graupel data
threshold: float
    threshold for graupel data
dilatation: int
    dilatation for graupel data
---------------------------------------------------------
OUT

---------------------------------------------------------
EXAMPLE CALL

---------------------------------------------------------
Killian P. Brennan
30.04.2024
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


def filter_hailcast_artefacts(hailcast, graupel, threshold=0.1, dilatation=40):
    # remove haicast data smaller than 10.1mm
    hailcast = hailcast.where(hailcast > 10.1, 0)
    hailcast_filtered = xr.Dataset()

    for t in range(0, len(hailcast.time) + 1):
        if t == 0:
            # get graupel data from hour after
            graupel_t = graupel.isel(time=slice(t, t + 2)).max(dim="time")
        elif t == len(hailcast.time):
            # get graupel data from hour before
            graupel_t = graupel.isel(time=slice(t - 1, t + 1)).max(dim="time")
        else:
            # get graupel data from hour before and after
            graupel_t = graupel.isel(time=slice(t - 1, t + 2)).max(dim="time")

        # apply threshold
        graupel_t = graupel_t.where(graupel_t > threshold, 0)

        # apply binary dilatation
        structure = disk(dilatation)
        graupel_t_dilated = ndimage.binary_dilation(graupel_t, structure=structure)

        # apply resulting mask
        hailcast_t = hailcast.isel(time=t)
        hailcast_t_filtered = hailcast_t.where(graupel_t_dilated, 0)

        hailcast_filtered = xr.concat([hailcast_filtered, hailcast_t_filtered], dim="time")

    return hailcast_filtered


def main(args):

    # load data
    hailcast = xr.open_dataset(args.hailcast)
    graupel = xr.open_dataset(args.graupel)

    # filter artefacts
    hailcast_filtered = filter_hailcast_artefacts(
        hailcast, graupel, args.threshold, args.dilatation
    )

    # save data
    hailcast_filtered.to_netcdf(args.output)


if __name__ == "__main__":

    parser = ArgumentParser(description="filter orographic hailcast artefacts")
    parser.add_argument("--hailcast", type=str, help="path to hailcast data")
    parser.add_argument("--graupel", type=str, help="path to graupel data")
    parser.add_argument("--threshold", type=float, help="threshold for hailcast data")
    parser.add_argument("--dilatation", type=int, help="dilatation for hailcast data")
    parser.add_argument("--output", type=str, help="path to output data")

    args = parser.parse_args()

    main(args)
