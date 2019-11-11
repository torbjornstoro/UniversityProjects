from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
from scipy.stats import norm
import pims
import trackpy as tp

# Tweaking matplotlib styles
mpl.rc('figure',  figsize=(10, 6))
mpl.rc('image', cmap='gray')
#plt.style.use("ggplot")

# Set path to umbrella result folder
gen_path = "/Users/petter/Desktop/TFY4335_LAB"

# Specify the relevant cases to iterate through
cases = {
    "A_BF": ["/Sample_A_BF_40x_frames/*.png",
             "/Sample_A_BF_40x_result_plots/",
             "A_BF_40x"],
    "B_BF": ["/Sample_B_BF_40x_frames/*.png",
             "/Sample_B_BF_40x_result_plots/",
             "B_BF_40x"],
    "A_DF": ["/Sample_A_DF_40x_frames/*.png",
             "/Sample_A_DF_40x_result_plots/",
             "A_DF_40x"],
    "B_DF": ["/Sample_B_DF_40x_frames/*.png",
             "/Sample_B_DF_40x_result_plots/",
             "B_DF_40x"]
}

# Set True/False to run respective code sections
plots = {
    "Generate_batch":               False,
    "Annotate_filtered":            True,
    "Annotate_unfiltered":          True,
    "Size_vs_mass":                 True,
    "Hydrodynamic_radius_filtered": True,
    "Gyration_radius_unfiltered":   False,
    "Gyration_radius_filtered":     False,
    "Trajectory_drift":             False,
    "Variance_linear":              True,
    "Variance_power_fit":           True,
    "Variance_all_parts":           True,
    "Return_A_and_D":               False
}

# Defining some relevant parameters
k_b         = 1.38*10**-23
T           = 298.15
μ           = 10**-3
ratio_μm_px = 308/1360.
fps         = 7.5
size        = 20

# Main function for running specified code blocks
def main():

    for case_idx, case in enumerate(cases.values()):

        res_path = gen_path + case[1]
        frames = pims.ImageSequence(gen_path + case[0], as_grey=True)

        # Stores the unfiltered annotated image of a frame to local file path
        if plots["Annotate_unfiltered"]:
            k = tp.locate(frames[0], 11, invert=[case_idx in [0, 1]], minmass=200)
            fig = plt.figure("Annotated_unfiltered_image_" + case[2])
            ax1 = fig.add_subplot()
            a = ["k", "w"][case_idx in [2, 3]]
            tp.annotate(k, frames[0], color=a, ax=ax1)
            #ax1.set_title("Annotated unfiltered image case: "+ case[2])
            #ax1.set_xlabel("x[px]", fontsize=size)
            #ax1.set_ylabel("y[px]", fontsize=size)
            ax1.tick_params(axis="both", which="both", top=False, bottom=False, labelbottom=False, right=False,
                            left=False, labelleft=False)
            fig.savefig(res_path + "Annotated_unfiltered_image_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close(fig)

        # If True: Tracks all frames and stores .csv to locally, else: imports such .csv file from local
        if plots["Generate_batch"]:
            f = tp.batch(frames[:225], 11, minmass=100, invert=[case_idx in [0, 1]])
            f.to_csv(res_path + "batch_" + case[2] + ".csv")
        if not plots["Generate_batch"]:
            f = pd.read_csv(res_path + "batch_" + case[2] + ".csv")

        # Linking and filtering
        t = tp.link_df(f, 5, memory=3)
        t1 = tp.filter_stubs(t, 50)

        # Plots the size vs mass profile and saves to local file path
        if plots["Size_vs_mass"]:
            fig = plt.figure("Size_vs_mass_" + case[2])
            ax1 = fig.add_subplot()
            tp.mass_size(t1.groupby('particle').mean(), ax=ax1)# convenience function -- just plots size vs. mass
            #ax1.set_title("Size vs mass case: " + case[2])
            ax1.set_xlabel("mass", fontsize=size)
            ax1.set_ylabel("Gyration radius [px]", fontsize=size)
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.tick_params(axis="both", which="major", labelsize=size)
            ax1.tick_params(axis="both", which="minor", labelsize=size)
            fig.savefig(res_path + "Size_vs_mass_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close(fig)

        if plots["Annotate_filtered"]:

            if case_idx in [0, 1]:  # Set BF condition
                condition = lambda x: ((x['mass'].mean() > 250) & (x['size'].mean() < 3.0) &
                                       (x['ecc'].mean() < 0.1))

            elif case_idx in [2, 3]:  # Set DF condition
                condition = lambda x: ((x['mass'].mean() > 100) & (x['size'].mean() < 5.0) &
                                       (x['ecc'].mean() < 0.1))

            t2 = tp.filter(t1, condition)  # a wrapper for pandas' filter that works around a bug in v 0.12

            fig = plt.figure("Annotated_filtered_image_" + case[2])
            ax1 = fig.add_subplot()
            k = ["k", "w"][case_idx in [2, 3]]
            tp.annotate(t2[t2['frame'] == 0], frames[0], color=k, ax=ax1)
            #ax1.set_title("Annotated filtered image case: " + case[2])
            #ax1.set_xlabel("x[px]", fontsize=size)
            #ax1.set_ylabel("y[px]", fontsize=size)
            ax1.tick_params(axis="both", which="both", top=False, bottom=False, labelbottom=False, right=False,
                            left=False, labelleft=False)
            fig.savefig(res_path + "Annotated_filtered_image_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close(fig)

        if plots["Gyration_radius_filtered"]:
            size_dis_t1 = [i * ratio_μm_px for i in t1['size']]
            plt.figure("Gyration_radius_filtered_" + case[2])
            plt.hist(size_dis_t1, bins=300, color="k", alpha=0.5)
            #plt.title("Gyration radius filtered case: "+ case[2])
            plt.ylabel("Events", fontsize=size)
            plt.xlabel("Gyration radius [μm]", fontsize=size)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.tick_params(axis="both", which="major", labelsize=size)
            plt.tick_params(axis="both", which="minor", labelsize=size)
            plt.savefig(res_path + "Gyration_radius_filtered" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close("all")

        if plots["Gyration_radius_unfiltered"]:
            size_dis_t = [i * ratio_μm_px for i in t['size']]
            plt.figure("Gyration_radius_unfiltered_" + case[2])
            plt.hist(size_dis_t, bins=300, color="k", alpha=0.5)
            #plt.title("Gyration radius unfiltered case: " + case[2])
            plt.ylabel("Events", fontsize=size)
            plt.xlabel("Gyration radius [μm]", fontsize=size)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.tick_params(axis="both", which="major", labelsize=size)
            plt.tick_params(axis="both", which="minor", labelsize=size)
            plt.savefig(res_path + "Gyration_radius_unfiltered" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close("all")

        d = tp.compute_drift(t1)
        tm = tp.subtract_drift(t1, d)

        if plots["Trajectory_drift"]:
            fig = plt.figure("Trajectory_drift_subtracted_" + case[2])
            ax1 = fig.add_subplot()
            ax1.tick_params(axis="both", which="major", labelsize=size)
            ax1.tick_params(axis="both", which="minor", labelsize=size)
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            tp.plot_traj(tm, ax=ax1)
            #ax1.set_title("Trajectory with drift subtracted case: " + case[2])
            plt.savefig(res_path + "Trajectory_drift_subtracted_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close(fig)

        if plots["Variance_all_parts"]:
            im = tp.imsd(tm, ratio_μm_px, fps, max_lagtime=225)  # microns per pixel = 100/285., frames per second = 24
            plt.figure("Variance_for_all_particles_" + case[2])
            #plt.title("Variance for all particles case: " + case[2])
            plt.plot(im.index, im, 'k-', alpha=0.1)  # black lines, semitransparent
            plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]', fontsize=size),
            plt.xlabel('time $t$', fontsize=size)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.xscale('log')
            plt.yscale('log')
            plt.savefig(res_path + "Variance_for_all_particles_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)

        if plots["Variance_linear"] or plots["Variance_power_fit"] or plots["Return_A_and_D"]:
            em = tp.emsd(tm, ratio_μm_px, fps, max_lagtime=225)

        if plots["Variance_linear"]:
            plt.figure("Variance_linear_fit_" + case[2])
            #plt.title("Variance linear fit case: " + case[2])
            plt.plot(em.index, em, 'ko', alpha=0.5)
            plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]', fontsize=size),
            plt.xlabel('time $t$ [s]', fontsize=size)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.tick_params(axis="both", which="major", labelsize=size)
            plt.tick_params(axis="both", which="minor", labelsize=size)
            plt.xscale('log')
            plt.yscale('log')
            plt.ylim(1e-2, 50)
            plt.savefig(res_path + "Variance_linear_fit_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)

        if plots["Variance_power_fit"]:
            fig = plt.figure("Variance_power_fit_" + case[2])
            ax1 = fig.add_subplot()
            tp.utils.fit_powerlaw(em, ax=ax1, color="k", alpha=0.5)
            #ax1.set_title("Variance power fitted case: " + case[2])
            ax1.set_ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]', fontsize=size)
            ax1.set_xlabel('time $t$ [s]', fontsize=size)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            ax1.tick_params(axis="both", which="major", labelsize=size)
            ax1.tick_params(axis="both", which="minor", labelsize=size)
            fig.savefig(res_path + "Variance_power_fit_" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close(fig)

        if plots["Hydrodynamic_radius_filtered"]:
            im = tp.imsd(tm, ratio_μm_px, fps, max_lagtime=225)
            r_h = []
            count = 0
            im = im.rename_axis("ID").values

            for index in range(1, len(im[0])):
                if isinstance(im[40][index], float) and isinstance(im[8][index], float):
                    D = (im[40][index] - im[8][index])/(4 * (40-8)/fps) * 10 **(-12)
                    if isinstance(D, float):
                        r_h += [abs(10**6*(k_b * T) / (6 * np.pi * μ * D))]
                        if 0 < abs(10**6*(k_b * T) / (6 * np.pi * μ * D)) < 6:
                            count += 1

            print("In interval: ", count, "Total: ", len(r_h), "Ratio: ", count/len(r_h))
            plt.figure("Hydrodynamic_radius_filtered_" + case[2])
            plt.hist(r_h, bins=int(count/3), color="k", alpha=0.5, range=(0, 6))
            #plt.title("Hydrodynamic radius filtered case: "+ case[2])
            plt.ylabel("Trajectories", fontsize=size)
            plt.xlabel("Hydrodynamic radius [μm]", fontsize=size)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.tick_params(axis="both", which="major", labelsize=size)
            plt.tick_params(axis="both", which="minor", labelsize=size)
            plt.savefig(res_path + "Hydrodynamic_radius_filtered" + case[2] + ".png", bbox_inches="tight", pad_inches=0)
            plt.close("all")

        if plots["Return_A_and_D"]:
            A = tp.utils.fit_powerlaw(em, ax=ax1, color="k", alpha=0.5)["A"]*10**(-12)
            print(tp.utils.fit_powerlaw(em, ax=ax1, color="k", alpha=0.5))
            print("For case ", case[2], " A=", A, ", and D=", A/4, ".")

if __name__ == "__main__":
    main()