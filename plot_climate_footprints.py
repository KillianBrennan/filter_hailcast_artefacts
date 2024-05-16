'''
---------------------------------------------------------
EXAMPLE CALL
python /home/kbrennan/filter_hailcast_artefacts/plot_climate_footprints.py /net/litho/atmosdyn2/kbrennan/data/climate/present/ /net/litho/atmosdyn2/kbrennan/data/climate/present/plots/potential_artefacts/ 20210102 20221231
---------------------------------------------------------
'''


import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd

from argparse import ArgumentParser

def main(inpath, outpath, start_day, end_day):

    # make sure the output directory exists
    os.makedirs(outpath, exist_ok=True)

    const_path = os.path.join(inpath, "lffd20101019000000c.nc")
    const = xr.open_dataset(const_path)

    origpath = os.path.join(inpath, "5min_2D")
    filteredpath = os.path.join(inpath, "filtered_hailcast")

    start_day = pd.to_datetime(start_day, format="%Y%m%d")
    end_day = pd.to_datetime(end_day, format="%Y%m%d")
    days = pd.date_range(start_day, end_day, freq="D")
    for day in days:
        day_str = day.strftime("%Y%m%d")
        print(day_str)

        filtered = xr.open_dataset(os.path.join(filteredpath, 'lffd' + day_str + "_0606.nc"))
        orig = xr.open_dataset(os.path.join(origpath, 'lffd' + day_str + "_0606.nz"))

        orig_footprint = orig.max(dim="time")
        filtered_footprint = filtered.max(dim="time")

        diff = orig_footprint.DHAIL_MX - filtered_footprint.DHAIL_MX

        # n_artefacts = np.sum(diff.values > 5)
        # plot
        fig, ax = plt.subplots(1, 1, figsize=(20, 16))
        pmesh = ax.pcolormesh(
            orig_footprint.rlon,
            orig_footprint.rlat,
            orig_footprint.TOT_PREC *12,
            vmax=20,
            vmin=0,
            cmap="Blues",
        )
        ax.contour(
            const.rlon,
            const.rlat,
            const.HSURF.squeeze(),
            levels=[1000,2000],
            colors=[(0, 0, 0, 0.5), (0, 0, 0, 1)],
            linewidths=0.5,
        )
        ax.contour(
            orig_footprint.rlon,
            orig_footprint.rlat,
            orig_footprint.DHAIL_MX,
            levels=[5],
            colors="red",
            linewidths=0.5,
        )

        # print number of pixels with artefacts
        # if n_artefacts > 0:
        #     print("n_artefacts: " + str(n_artefacts))
        # # add text to plot
        # ax.text(
        #     0.01,
        #     0.99,
        #     "n_artefacts: " + str(n_artefacts),
        #     transform=ax.transAxes,
        #     fontsize=10,
        #     verticalalignment="top",
        #     horizontalalignment="left",
        # )

        ax.contourf(
            diff.rlon,
            diff.rlat,
            diff,
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
        plt.savefig(os.path.join(outpath, day_str + ".png"), dpi=600)


if __name__ == "__main__":

    parser = ArgumentParser(description="filter orographic hailcast artefacts")

    parser.add_argument("inpath", type=str, help="path to the input directory")
    parser.add_argument("outpath", type=str, help="path to the output directory")
    parser.add_argument("start_day", type=str, help="start day in format YYYYMMDD")
    parser.add_argument("end_day", type=str, help="end day in format YYYYMMDD")

    args = parser.parse_args()

    main(args.inpath, args.outpath, args.start_day, args.end_day)



# import os

# import matplotlib.pyplot as plt
# import numpy as np
# import xarray as xr


# def main():
#     root_path = "/home/kbrennan/phd/data/climate/present"
#     footprint_path = os.path.join(root_path, "hail_footprints")
#     horizontal_path = os.path.join(root_path, "5min_2D")
#     filtered_path = os.path.join(root_path, "filtered_hailcast")
#     const_path = os.path.join(root_path, "lffd20101019000000c.nc")
#     const = xr.open_dataset(const_path)

#     footprint_files = os.listdir(footprint_path)
#     footprint_files.sort()
#     # only nc files
#     footprint_files = [f for f in footprint_files if f.endswith(".nc")]
#     # only the ones that start with ffd2021
#     footprint_files = [f for f in footprint_files if f.startswith("lffd2021")]
#     footprint_files = footprint_files[11:]
#     plotdir = os.path.join(root_path, "plots", "potential_artefacts")
#     if not os.path.exists(plotdir):
#         os.makedirs(plotdir)

#     for f in footprint_files:
#         day = f.split("_")[0][4:]
#         print(day)
#         footprint = xr.open_dataset(os.path.join(footprint_path, f))
#         filtered = xr.open_dataset(
#             os.path.join(filtered_path, "lffd" + day + "_0606.nc")
#         )

#         try:
#             horizontal = xr.open_mfdataset(
#                 os.path.join(horizontal_path, "*" + day + "*.nz")
#             )
#         except:
#             print("No file found for " + day)
#             continue

#         diff = footprint.DHAIL_MX.squeeze() - filtered.DHAIL_MX_filtered.max(dim="time")
#         n_artefacts = np.sum(diff.values > 5)
#         # plot
#         fig, ax = plt.subplots(1, 1, figsize=(20, 16))
#         pmesh = ax.pcolormesh(
#             horizontal.rlon,
#             horizontal.rlat,
#             horizontal.TOT_PREC.max(dim="time") * 60 / 5,
#             vmax=20,
#             vmin=0,
#             cmap="Blues",
#         )
#         ax.contour(
#             const.rlon,
#             const.rlat,
#             const.HSURF.squeeze(),
#             levels=[2000],
#             colors="black",
#             linewidths=0.5,
#         )
#         ax.contour(
#             footprint.rlon,
#             footprint.rlat,
#             footprint.DHAIL_MX.squeeze(),
#             levels=[5],
#             colors="red",
#             linewidths=0.5,
#         )

#         # print number of pixels with artefacts
#         if n_artefacts > 0:
#             print("n_artefacts: " + str(n_artefacts))
#         # add text to plot
#         ax.text(
#             0.01,
#             0.99,
#             "n_artefacts: " + str(n_artefacts),
#             transform=ax.transAxes,
#             fontsize=10,
#             verticalalignment="top",
#             horizontalalignment="left",
#         )

#         ax.contourf(
#             filtered.rlon,
#             filtered.rlat,
#             diff,
#             levels=[5, 50],
#             colors="tab:orange",
#             alpha=0.5,
#         )

#         # add colorbar
#         cbar = fig.colorbar(pmesh, ax=ax)
#         cbar.set_label("max precipitation rate (mm/h)")

#         # limit to swiss lat /lon
#         # ax.set_xlim(5, 11)
#         # ax.set_ylim(45, 48)

#         ax.set_title(day)
#         plt.savefig(os.path.join(plotdir, day + ".png"), dpi=600)


# if __name__ == "__main__":
#     main()
