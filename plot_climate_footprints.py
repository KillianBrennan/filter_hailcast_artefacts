import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def main():
    root_path = "/home/kbrennan/phd/data/climate/present"
    footprint_path = os.path.join(root_path, "hail_footprints")
    horizontal_path = os.path.join(root_path, "1h_2D")
    const_path = os.path.join(root_path, "lffd20101019000000c.nc")
    const = xr.open_dataset(const_path)

    footprint_files = os.listdir(footprint_path)
    footprint_files.sort()
    # only nc files
    footprint_files = [f for f in footprint_files if f.endswith(".nc")]
    # only the ones that start with ffd2021
    footprint_files = [f for f in footprint_files if f.startswith("lffd2021")]
    plotdir = os.path.join(root_path, "plots", "potential_artefacts")
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    for f in footprint_files:
        day = f.split("_")[0][4:]
        print(day)
        footprint = xr.open_dataset(os.path.join(footprint_path, f))

        horizontal = xr.open_mfdataset(
            os.path.join(horizontal_path, "*" + day + "*.nc")
        )

        # plot
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        pmesh = ax.pcolormesh(
            horizontal.lon,
            horizontal.lat,
            horizontal.TQG.max(dim="time"),
            vmax=5,
            vmin=0,
            cmap="Purples",
        )
        ax.contour(
            const.lon,
            const.lat,
            const.HSURF.squeeze(),
            levels=[2000, 3000],
            colors="black",
            linestyles=["dashed", "solid"],
        )
        ax.contour(
            footprint.lon,
            footprint.lat,
            footprint.DHAIL_MX.squeeze(),
            levels=[1],
            colors="red",
        )
        # add colorbar
        cbar = fig.colorbar(pmesh, ax=ax)
        cbar.set_label("TQG [kg/m^2]")

        # limit to swiss lat /lon
        ax.set_xlim(5, 11)
        ax.set_ylim(45, 48)

        ax.set_title(day)
        plt.savefig(os.path.join(plotdir, day + ".png"))


if __name__ == "__main__":
    main()
