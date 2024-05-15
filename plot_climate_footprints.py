import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def main():
    root_path = "/home/kbrennan/phd/data/climate/present"
    footprint_path = os.path.join(root_path, "hail_footprints")
    horizontal_path = os.path.join(root_path, "5min_2D")
    filtered_path = os.path.join(root_path, "filtered_hailcast")
    const_path = os.path.join(root_path, "lffd20101019000000c.nc")
    const = xr.open_dataset(const_path)

    footprint_files = os.listdir(footprint_path)
    footprint_files.sort()
    # only nc files
    footprint_files = [f for f in footprint_files if f.endswith(".nc")]
    # only the ones that start with ffd2021
    footprint_files = [f for f in footprint_files if f.startswith("lffd2021")]
    # footprint_files = footprint_files[-20:]
    plotdir = os.path.join(root_path, "plots", "potential_artefacts")
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    for f in footprint_files:
        day = f.split("_")[0][4:]
        print(day)
        footprint = xr.open_dataset(os.path.join(footprint_path, f))
        filtered = xr.open_dataset(
            os.path.join(filtered_path, "lffd" + day + "_0606.nc")
        )

        try:
            horizontal = xr.open_mfdataset(
                os.path.join(horizontal_path, "*" + day + "*.nz")
            )
        except:
            print("No file found for " + day)
            continue

        # plot
        fig, ax = plt.subplots(1, 1, figsize=(20, 16))
        pmesh = ax.pcolormesh(
            horizontal.rlon,
            horizontal.rlat,
            horizontal.TOT_PREC.max(dim="time") * 60 / 5,
            vmax=20,
            vmin=0,
            cmap="Blues",
        )
        ax.contour(
            const.rlon,
            const.rlat,
            const.HSURF.squeeze(),
            levels=[2000],
            colors="black",
            linewidths=0.5,
        )
        ax.contour(
            footprint.rlon,
            footprint.rlat,
            footprint.DHAIL_MX.squeeze(),
            levels=[1],
            colors="red",
            linewidths=0.5,
        )

        ax.contourf(
            filtered.rlon,
            filtered.rlat,
            filtered.DHAIL_MX.max(dim="time") - footprint.DHAIL_MX.squeeze(),
            levels=[5, 50],
            colors="tab:orange",
            alpha=0.5,
        )
        # add colorbar
        cbar = fig.colorbar(pmesh, ax=ax)
        cbar.set_label("max precipitation rate (mm/h)")

        # limit to swiss lat /lon
        # ax.set_xlim(5, 11)
        # ax.set_ylim(45, 48)

        ax.set_title(day)
        plt.savefig(os.path.join(plotdir, day + ".png"), dpi=600)


if __name__ == "__main__":
    main()
