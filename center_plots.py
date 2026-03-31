import numpy as np
import matplotlib.ticker as tick
import matplotlib.pyplot as plt

from utilities import raFormatter, decFormatter
from center_methods import getCenters
from astro_data_format import toAstroData

# Plottar varje graf från centrum och formatterar axlar
def plotter(ax, centers, color, label):
    ax.plot(centers[1], centers[2], marker='*', color=color, markerfacecolor="gold", markeredgecolor="darkorange", label=label, markersize=3.0, alpha=0.5)
    ax.axis('equal') # 1 grad RA = 1 grad dec
    ax.set_xlabel("Ra")
    ax.set_ylabel("Dec")

    # Formattera ticks så att de har rimligt avstånd (för den långa texten) och rätt dec/ra format
    max_scale = max(np.max(centers[1]) - np.min(centers[1]), np.max(centers[2]) - np.min(centers[2]))
    ax.xaxis.set_major_locator(tick.MultipleLocator(base=0.4*max_scale))
    ax.yaxis.set_major_locator(tick.MultipleLocator(base=0.2*max_scale))
    ax.xaxis.set_major_formatter(tick.FuncFormatter(raFormatter))
    ax.yaxis.set_major_formatter(tick.FuncFormatter(decFormatter))


def plotMethods(method1, method1_str, method2, method2_str, R_Dor, R_Leo, W_Hya):

    # Beräkna centrum för de två metoderna för varje stjärna
    print("R Dorados", method1_str, "\n=============================")
    R_Dor_center1 = getCenters(method1, R_Dor)
    print("R Dorados", method2_str, "\n=============================")
    R_Dor_center2 = getCenters(method2, R_Dor)
    print("R Leonis", method1_str, "\n=============================")
    R_Leo_center1 = getCenters(method1, R_Leo)
    print("R Leonis", method2_str, "\n=============================")
    R_Leo_center2 = getCenters(method2, R_Leo)
    print("W Hydrae", method1_str, "\n=============================")
    W_Hya_center1 = getCenters(method1, W_Hya)
    print("W Hydrae", method2_str, "\n=============================")
    W_Hya_center2 = getCenters(method2, W_Hya)

    # Spara data:
    toAstroData(fname = "R_Dor_1.astrom.dat", header=f"name = R Dorados - {method1_str}", data=R_Dor_center1)
    toAstroData(fname = "R_Dor_2.astrom.dat", header=f"name = R Dorados - {method2_str}", data=R_Dor_center2)
    toAstroData(fname = "R_Leo_1.astrom.dat", header=f"name = R Leonis - {method1_str}", data=R_Leo_center1)
    toAstroData(fname = "R_Leo_2.astrom.dat", header=f"name = R Leonis - {method2_str}", data=R_Leo_center2)
    toAstroData(fname = "W_Hya_1.astrom.dat", header=f"name = W Hydrae - {method1_str}", data=W_Hya_center1)
    toAstroData(fname = "W_Hya_2.astrom.dat", header=f"name = W Hydrae - {method2_str}", data=W_Hya_center2)

    print("Star plot", "\n=============================")
    # Skapa figur för att plotta de olika fallen
    fig, axs = plt.subplots(2, 2, figsize=(12, 12), dpi=120)

    # Plotta alla centrum
    axs[0,0].set_title("R Dorados")
    axs[1,0].set_title("R Leonis")
    axs[0,1].set_title("W Hydrae")

    plotter(axs[0,0], R_Dor_center1, "cornflowerblue", method1_str)
    plotter(axs[1,0], R_Leo_center1, "crimson", method1_str)
    plotter(axs[0,1], W_Hya_center1, "forestgreen", method1_str)

    plotter(axs[0,0], R_Dor_center2, "deepskyblue", method2_str)
    plotter(axs[1,0], R_Leo_center2, "indianred", method2_str)
    plotter(axs[0,1], W_Hya_center2, "lawngreen", method2_str)

    axs[0,0].legend()
    axs[1,0].legend()
    axs[0,1].legend()

    axs[1,1].set_axis_off()

    return R_Dor_center1, R_Dor_center2, R_Leo_center1, R_Leo_center2, W_Hya_center1, W_Hya_center2
