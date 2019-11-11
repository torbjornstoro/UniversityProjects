import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

font_size    = 14
filter_value = 0.001

gen_path = "/Users/petter/Desktop/TFY4335_LAB"
in_path  = "/Observed_radii_in"
res_path = "/Observed_radii_res"
spec_paths = {"0": "/BF_40x_Observed_SIze_Distribution_A",
              "1": "/BF_40x_Observed_SIze_Distribution_B",
              "2": "/DF_40x_Observed_Size_Distribution_A",
              "3": "/DF_40x_Observed_Size_Distribution_B"}

for i in range(len(spec_paths)):
    df = pd.read_csv(gen_path + in_path + spec_paths[str(i)] + ".csv")
    d_set = []
    for sub_arr in df.values:
        if np.sqrt(sub_arr[1]/np.pi) >= filter_value:
            d_set += [np.sqrt(sub_arr[1]/np.pi)]

    plt.figure()
    plt.hist(d_set, bins=30, color="k", alpha=0.5, rwidth=1)
    plt.tick_params(axis="both", which="both", labelsize=font_size)
    plt.ylabel("Trajectories", fontsize=font_size)
    plt.xlabel("Observed radius [μm]", fontsize=font_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(gen_path + res_path + spec_paths[str(i)] + ".png", bbox_inches="tight", pad_inches=0)

    print("The mean is: " + str(round(np.mean(d_set), 2)))