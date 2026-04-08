import numpy as np
import matplotlib.ticker as tick
import matplotlib.pyplot as plt

from utilities import raFormatter, decFormatter
from center_methods import getCenters
from astro_data_format import toAstroData, fromAstroData

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

# Det behövs färger om man ska plotta flera metoder
färger=["cornflowerblue", "crimson","forestgreen","deepskyblue","indianred","lawngreen", "aqua", "aquamarine", "black", "blueviolet", "chartreuse", "cyan" "darkgreen", "darkmagenta", "darksalmon", "deeppink", "goldenrod","hotpink", "indigo" , "lime", "olive", "seagreen"]
#Lista alla metoder och korresponderande filnamn
metoder={"maxIntensitet": ["R_Dor_MaxInt.astrom.dat","R_Leo_MaxInt.astrom.dat","W_Hya_MaxInt.astrom.dat"],
         "gaussIntensitet":["R_Dor_Gauss.astrom.dat","R_Leo_Gauss.astrom.dat","W_Hya_Gauss.astrom.dat",],
         "halfMax":["R_Dor_HalfMax.astrom.dat","R_Leo_HalfMax.astrom.dat","W_Hya_HalfMax.astrom.dat"],
         "halfViktad":["R_Dor_HalfViktad.astrom.dat","R_Leo_HalfViktad.astrom.dat","W_Hya_HalfViktad.astrom.dat"],
         "halfLSQ":["R_Dor_HalfLSQ.astrom.dat","R_Leo_HalfLSQ.astrom.dat","W_Hya_HalfLSQ.astrom.dat"],
         "moffat":["R_Dor_Moffat.astrom.dat","R_Leo_Moffat.astrom.dat","W_Hya_Moffat.astrom.dat"]
}

def plotMethods(method):
    # Skapa listor
    R_Dor_center=[]
    R_Leo_center=[]
    W_Hya_center=[]
    print("Star plot", "\n=============================")
    # Skapa figur för att plotta de olika fallen
    fig, axs = plt.subplots(2, 2, figsize=(12, 12), dpi=120)
    # Ge plottsen namn. 
    axs[0,0].set_title("R Dorados")
    axs[1,0].set_title("R Leonis")
    axs[0,1].set_title("W Hydrae") 
    # Kör igenom alla olika metoder
    for i in range(0, len(method)):
        # Hämta data för specifika metoder för varje stjärna
        R_Dor_center.append(fromAstroData(metoder.get(method[i])[0]))
        R_Leo_center.append(fromAstroData(metoder.get(method[i])[1]))
        W_Hya_center.append(fromAstroData(metoder.get(method[i])[2]))

        # Vi måste plotta böset
        plotter(axs[0,0], R_Dor_center[i], färger[i], method[i])
        plotter(axs[1,0], R_Leo_center[i], färger[i], method[i])
        plotter(axs[0,1], W_Hya_center[i], färger[i], method[i])
        axs[0,0].legend()
        axs[1,0].legend()
        axs[0,1].legend()

        axs[1,1].set_axis_off()
    return R_Dor_center,R_Leo_center, W_Hya_center
