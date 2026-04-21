import numpy as np
import scipy as sci
from astropy.modeling import models, fitting



# Temporary
import matplotlib.pyplot as plt

##########################
# MAIN Center function
##########################

# Få centerpunkter med metod och observationer
# Ger tillbaka tidpunkter, ra:s och dec:s
def getCenters(method, observations, lowLim=0.6, upLim=1, debug=False, halva=False, user=False):
    ts = []
    ras = []
    decs = []
    ras_err = []
    decs_err = []

    # Gå igenom varje observaton
    for observation in observations:
        ts.append(observation[1])
        if debug:
            print("Band: ", observation[0], "Tidpunkt: ", observation[1])
        if halva:
            ra, dec, intensitet, sigma_ra, sigma_dec = method(observation[2], observation[3], observation[4], lowLim, upLim, debug=debug)
        elif user:
            ra, dec, sigma_ra, sigma_dec = observation[7], observation[8], observation[5], observation[6]
        else:
            ra, dec, intensitet, sigma_ra, sigma_dec = method(observation[2], observation[3], observation[4], debug=debug)

        ras.append(ra)
        decs.append(dec)
        ras_err.append(sigma_ra)
        decs_err.append(sigma_dec)

    # Sortera datan enligt observationstid
    idx = np.argsort(ts)

    ts = np.array(ts)[idx]
    ras = np.array(ras)[idx]
    decs = np.array(decs)[idx]
    ras_err = np.array(ras_err)[idx]
    decs_err = np.array(decs_err)[idx]
    data = np.array([ts, ras, decs, ras_err, decs_err])

    return data

######################
# CENTER FUNCTIONS
######################

# Hitta centrum med ljusaste punkt:
def maxIntensitet(ras, decs, intensities, debug=False):
    # Hämta index för ljusaste intensitet
    j, i = np.unravel_index(np.argmax(intensities, axis=None), intensities.shape)

    return ras[i], decs[j], intensities[j, i], np.abs(ras[0]-ras[1]), np.abs(decs[0]-decs[1])

# Hitta centrum med gauss anpassning:
def gaussIntensitet(ras, decs, intensities, debug=False):
    r, d = np.meshgrid(ras, decs)

    # Få ut maximala intensiteten som utgångspunkt för den gaussiansk anpassningen
    ra_max, dec_max, int_max, _, _ = maxIntensitet(ras, decs, intensities)

    # Standard avvikelse satt till 12 arcsec, rimlig approximation av stjärnans radie
    stddev_factor = 12/(360*60*60)

    # Initala gauss modellen: Startpunkt för anpassningen
    gauss_init = models.Gaussian2D(amplitude=1, x_mean=ra_max, y_mean=dec_max, \
                                   x_stddev=stddev_factor, y_stddev=stddev_factor)
    
    # Anpassningsmetoden som fungerade bäst. TODO: Undersök varför och om det finns bättre anpassare.
    gauss_fitter = fitting.LevMarLSQFitter()

    # Gör gaussisk anpassning med en maximal iteration på 10'000, (lägre gav sämre resultat)
    gauss_model = gauss_fitter(gauss_init, r, d, intensities, maxiter=10000) 

    # Extrahera fel från metod
    residuals = intensities - gauss_model(r, d)
    N = intensities.size
    p = len(gauss_model.parameters)

    sigma2 = np.sum(residuals**2) / (N - p)
    cov = gauss_fitter.fit_info['param_cov']
    cov = cov * sigma2
    sigma_ra = np.sqrt(cov[1, 1])
    sigma_dec = np.sqrt(cov[2, 2]) 

    # För debug så plottas den initiala gauss modellen, den gaussiska anpassningen och ursprungliga stjärnbilden bredvid varandra.
    if debug:
        fig, axs = plt.subplots(1, 3, figsize=(12, 4), dpi=120)
        fig.suptitle("")
        axs[0].set_title("Init")
        axs[0].pcolormesh(ras, decs, gauss_init(r, d))
        axs[0].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[1].set_title("Model")
        axs[1].pcolormesh(ras, decs, gauss_model(r, d))
        axs[1].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[2].set_title("Actual")
        axs[2].pcolormesh(ras, decs, intensities)
        axs[2].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        plt.show()

    return gauss_model.x_mean.value, gauss_model.y_mean.value, gauss_model.amplitude.value, sigma_ra, sigma_dec

#Hittar ras och dec för ett intensitetsintervall
def hittaIntensitet(ras, decs, intensities, lowLim, upLim):

    #skaffa maxintensitet
    ra_max, dec_max, int_max, _, _ = maxIntensitet(ras, decs, intensities, debug=False)

    #skapa lite listor
    delint=[]
    delras=[]
    deldec=[]

    # Hitta data mellan två intensiteter
    for i in range(intensities.shape[0]):
        for j in range(intensities.shape[1]):
            if intensities[i,j]>lowLim*int_max and intensities[i,j]<upLim*int_max:
                delint.append(intensities[i,j])
                delras.append(ras[j])
                deldec.append(decs[i])

    return delint, delras, deldec

# Hitta centrum med halva max
def halfMax(ras, decs, intensities, lowLim=0.6, upLim=1, debug=False):

    #skaffa ras och dec för intensiteter i ett godtyckligt intervall mellan 0 och 1.
    delint, delras, deldec=hittaIntensitet(ras, decs, intensities, lowLim, upLim)

    # Medelvärdera, kan testa andra sätt att interpolera ett centrum
    meanRas=np.mean(delras)
    meanDec=np.mean(deldec)
    meanint=np.mean(delint)

    # Extrahera osäkerheter
    sigma_ra = np.std(delras, ddof = 1)/np.sqrt(len(delras))
    sigma_dec = np.std(deldec, ddof = 1)/np.sqrt(len(deldec))

    # Om man vill se figurer
    if debug:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4), dpi=120)
        fig.suptitle("")
        axs[0].set_title("Medelvärde")
        axs[0].pcolormesh(ras, decs, intensities)
        axs[0].plot(meanRas, meanDec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[1].set_title("Intensitetsintervall")
        axs[1].pcolormesh(ras, decs, intensities)
        axs[1].plot(delras, deldec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        plt.show()
    return meanRas, meanDec, meanint, sigma_ra, sigma_dec

# Hitta centrum med halva max viktad med intensiteter
def halfViktad(ras, decs, intensities, lowLim=0.6, upLim=1, debug=False):
    #skaffa ras och dec för intensiteter i ett godtyckligt intervall mellan 0 och 1.
    delint, delras, deldec=hittaIntensitet(ras, decs, intensities, lowLim, upLim)

    #medelvärdera med intensiteter som vikt
    meanRas=np.average(delras,weights=delint)
    meanDec=np.average(deldec, weights=delint)
    meanint=np.average(delint)

    #xtrahera osäkerheter
    w = np.array(delint)
    w_sum = np.sum(w)

    sigma_ra = np.sqrt(np.sum(w * (delras - meanRas)**2)) / w_sum
    sigma_dec = np.sqrt(np.sum(w * (deldec - meanDec)**2)) / w_sum
    # Om man vill se figurer
    if debug:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4), dpi=120)
        fig.suptitle("")
        axs[0].set_title("Viktat medelvärde")
        axs[0].pcolormesh(ras, decs, intensities)
        axs[0].plot(meanRas, meanDec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[1].set_title("Intensitetsintervall")
        axs[1].pcolormesh(ras, decs, intensities)
        axs[1].plot(delras, deldec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        plt.show()
        

    return meanRas, meanDec, meanint, sigma_ra, sigma_dec


# Hitta centrum med halva max LSQ, 
#R_Dor bäst vid 0.4-0.6
#R_Leo bäst vid 0.7-0.8
#W Hya bäst vid 0.5-0.7
def halfLSQ(ras, decs, intensities, lowLim=0.4, upLim=0.6, debug=False):
    #skaffa maxintensitet
    ra_max, dec_max, int_max, _, _ = maxIntensitet(ras, decs, intensities, debug=False)

    #skaffa ras och dec för intensiteter i ett godtyckligt intervall mellan 0 och 1.
    delint, delras, deldec=hittaIntensitet(ras, decs, intensities, lowLim, upLim)

    #Använd maxintensiteten som utgångspunkt 
    x_0=[ra_max,dec_max]

    #Beräkna kvadraterna på avstånden mellan x_0 och de valda intensiteterna
    def square(x_0):
        summa=sum((np.array(deldec)-x_0[1])**2+(np.array(delras)-x_0[0])**2)
        return summa

    #minimera square
    x_1=sci.optimize.minimize(square,x_0)
    
    meanint=np.average(delint)
    lsqRas=x_1.get("x")[0]
    lsqDec=x_1.get("x")[1]

    sigma_ra = np.std(delras, ddof = 1)/np.sqrt(len(delras))
    sigma_dec = np.std(deldec, ddof = 1)/np.sqrt(len(deldec))
    # Om man vill se figurer
    if debug:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4), dpi=120)
        fig.suptitle("")
        axs[0].set_title("Minsta kvadrat centrum")
        axs[0].pcolormesh(ras, decs, intensities)
        axs[0].plot(lsqRas, lsqDec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[1].set_title("Intensitetsintervall")
        axs[1].pcolormesh(ras, decs, intensities)
        axs[1].plot(delras, deldec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        plt.show()
        

    return lsqRas, lsqDec, meanint, sigma_ra, sigma_dec


# Hitta centrum med moffat anpassning 
def moffat(ras, decs, intensities, debug=False):
    r, d = np.meshgrid(ras, decs)
    #utgå från maxintensitet
    ra_max, dec_max, int_max, _, _ = maxIntensitet(ras, decs, intensities)

    stddev_factor = 12/(360*60*60)
    #modell moffat
    moffat_init = models.Moffat2D(amplitude=1, x_0=ra_max, y_0=dec_max,
                                   gamma=stddev_factor*2, alpha=2)
    #Samma fitter som fungerade bra för gauss
    moffat_fitter = fitting.LevMarLSQFitter()

    moffat_model = moffat_fitter(moffat_init, r, d, intensities, maxiter=1000)

    # Extrahera fel
    residuals = intensities - moffat_model(r, d)
    N = intensities.size
    p = len(moffat_model.parameters)

    sigma2 = np.sum(residuals**2) / (N - p)
    cov = moffat_fitter.fit_info['param_cov']
    cov = cov * sigma2
    sigma_ra = np.sqrt(cov[1, 1])
    sigma_dec = np.sqrt(cov[2, 2])  
    #om du vill se figurer
    if debug:
        fig, axs = plt.subplots(1, 3, figsize=(12, 4), dpi=120)
        fig.suptitle("")
        axs[0].set_title("Init")
        axs[0].pcolormesh(ras, decs, moffat_init(r, d))
        axs[0].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[1].set_title("Model")
        axs[1].pcolormesh(ras, decs, moffat_model(r, d))
        axs[1].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[2].set_title("Actual")
        axs[2].pcolormesh(ras, decs, intensities)
        axs[2].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        plt.show()
    
    return moffat_model.x_0.value, moffat_model.y_0.value, moffat_model.amplitude.value, sigma_ra, sigma_dec


# Något knäppt med denna, helt wack just nu
def IntIntegrering(ras, decs, intensities, debug=False):
    r, d = np.meshgrid(ras, decs)

    #Filtrera bort bakgrundsbrus genom att endast ta ett litet område runt den ljusaste punkten
    iy, ix = np.unravel_index(np.argmax(intensities), intensities.shape)

    box = 10

    sub = intensities[iy-box:iy+box+1, ix-box:ix+box+1]

    d_sub = decs[iy-box:iy+box+1]
    r_sub = ras[ix-box:ix+box+1]

    r_grid, d_grid = np.meshgrid(r_sub, d_sub)
    # Beräkna tyngdpunkten av intensiteterna
    weights = sub
    I_sum = np.sum(weights)
    x_iwc = np.sum(r_grid*weights) / I_sum
    y_iwc = np.sum(d_grid*weights) / I_sum

    # Beräkna osäkerheter
    sigma_ra = np.sqrt(np.sum(weights * (r_grid - x_iwc)**2) / I_sum)
    sigma_dec = np.sqrt(np.sum(weights * (d_grid - y_iwc)**2) / I_sum)
    N_eff = (I_sum**2) / np.sum(weights**2)
    sigma_ra /= np.sqrt(N_eff)
    sigma_dec /= np.sqrt(N_eff)

    # debug plots
    if debug:
        fig, axs = plt.subplots(figsize=(12, 4), dpi=120)
        fig.suptitle("")
        plt.title("Intensitetsintegrering centrum")
        plt.pcolormesh(ras, decs, intensities)
        plt.errorbar(x_iwc, y_iwc, xerr=sigma_ra, yerr=sigma_dec, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        print(sigma_ra, sigma_dec)
        plt.show()
    return x_iwc, y_iwc, intensities, sigma_ra, sigma_dec


# Hitta centrum med Lorentz2D anpassning
def Lorentz(ras, decs, intensities, debug=False):
    r, d = np.meshgrid(ras, decs)
    #Utgår från maxintensitet
    ra_max, dec_max, int_max, _, _ = maxIntensitet(ras, decs, intensities)

    #Lorentz2D model
    Lorentz_init = models.Lorentz2D(amplitude=int_max, x_0=ra_max, y_0=dec_max, fwhm=12/(360*60*60))

    Lorentz_fitter = fitting.DogBoxLSQFitter()

    Lorentz_model = Lorentz_fitter(Lorentz_init, r, d, intensities, maxiter=1000)

    # Extrahera fel
    residuals = intensities - Lorentz_model(r, d)
    N = intensities.size
    p = len(Lorentz_model.parameters)

    sigma2 = np.sum(residuals**2) / (N - p)
    cov = Lorentz_fitter.fit_info['param_cov']
    cov = cov * sigma2
    sigma_ra = np.sqrt(cov[1, 1])
    sigma_dec = np.sqrt(cov[2, 2])  

    #om du vill se figurer (Kopia av moffat)
    if debug:
        fig, axs = plt.subplots(1, 4, figsize=(16, 4), dpi=120)
        fig.suptitle("")
        axs[0].set_title("Init")
        axs[0].pcolormesh(ras, decs, Lorentz_init(r, d))
        axs[0].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[1].set_title("Model")
        axs[1].pcolormesh(ras, decs, Lorentz_model(r, d))
        axs[1].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[2].set_title("Actual")
        axs[2].pcolormesh(ras, decs, intensities)
        axs[2].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        axs[3].set_title("Residual")
        axs[3].pcolormesh(ras, decs, Lorentz_model(r,d)- intensities )
        axs[3].plot(ra_max, dec_max, marker='*', markerfacecolor="gold", markeredgecolor="darkorange", alpha=0.5, markersize=3.0)
        plt.show()

    return Lorentz_model.x_0.value, Lorentz_model.y_0.value, Lorentz_model.amplitude.value , sigma_ra, sigma_dec

