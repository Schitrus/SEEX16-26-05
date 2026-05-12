import numpy as np
import scipy as sci
import astropy as aspy
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
import jplephem # Need to install
from astropy.coordinates import solar_system_ephemeris
import matplotlib.gridspec as gridspec


solar_system_ephemeris.set("DE440")

# For calibration stars
def read_calibration_star(file):
    data_full = pd.read_csv(file, sep = r'\s+', header=None, names=['time', 'ra', 'ra_err', 'dec', 'dec_err'], usecols=[0,1,2,3,4])
    data = data_full.iloc[1:]
    refpos = data_full.iloc[0]
    return data, refpos

# Use astropy to get earths barycentric coordinates at time t
def bary_coords(t, ref_t):

    date = aspy.time.Time(t + ref_t, format = 'decimalyear')
    Epos = aspy.coordinates.get_body_barycentric('earth', date)

    X = Epos.x.to(u.au).value
    Y = Epos.y.to(u.au).value
    Z = Epos.z.to(u.au).value
    return X, Y, Z

# our model of the motion
def model(t, ref_t, *params):

    ra0, dec0, parallax, pm_ra, pm_dec = params
    X, Y, Z = bary_coords(t, ref_t)

    ra_model = ra0 + (pm_ra/np.cos(dec0))*t + (parallax/np.cos(dec0))*(X*np.sin(ra0) - Y*np.cos(ra0))
    dec_model = dec0 + pm_dec*t + parallax*(X*np.cos(ra0)*np.sin(dec0) + Y*np.sin(ra0)*np.sin(dec0) - Z*np.cos(dec0))

    return ra_model, dec_model

# Same as model, but without proper motion
def model_wopm(t, ref_t, *params):

    ra0, dec0, parallax, pm_ra, pm_dec = params
    
    X, Y, Z = bary_coords(t, ref_t)

    ra_model = ra0 + (parallax/np.cos(dec0))*(X*np.sin(ra0) - Y*np.cos(ra0))
    dec_model = dec0 + parallax*(X*np.cos(ra0)*np.sin(dec0) + Y*np.sin(ra0)*np.sin(dec0) - Z*np.cos(dec0))

    return ra_model, dec_model

# Take observations, and remove its proper motion
def model_subtract_pm(t, ras, decs, *params):
    
    ra0, dec0, parallax, pm_ra, pm_dec = params
        
    ra_sub_pm = ras - ((pm_ra/np.cos(dec0)) * t)
    dec_sub_pm = decs - (pm_dec * t)
    #print(f'ras before: {ras}, ras after {ra_sub_pm}')
    #print(f'decs before: {decs}, decs after {dec_sub_pm}')
    return ra_sub_pm, dec_sub_pm

# residuals to use with least square fitting
def residuals(params, t, ref_t, ras, decs, ras_err, decs_err):
    ra_model, dec_model = model(t, ref_t, *params)
    dec0 = params[1]
    dra = np.arctan2(
    np.sin(ras - ra_model),
    np.cos(ras - ra_model)
)
    residual_ra = dra / ras_err
    residual_dec = (decs - dec_model) / decs_err

    return np.concatenate([residual_ra, residual_dec])

# Function to fit for parameters
def fit_model(t, ref_t, ras, decs, initial, bounds, ras_err, decs_err):

    result = sci.optimize.least_squares(residuals, initial,  args = (t, ref_t, ras, decs, ras_err, decs_err), max_nfev=50000, bounds = bounds)

    ndata = 2 * len(ras)
    nparams = len(result.x)

    # Chi-square
    chi2 = 2 * result.cost

    # Reduced chi-square
    chi2_red = chi2 / (ndata - nparams)

    print(f"Reduced X^2 = {chi2_red:.3f}")
    return result

# Function for plotting
def makeplots(name, result, t, ref_t, ras, decs, ras_err, decs_err, save = False, datum=False):

    # Handle uncertainties
    ras_err = 3.6e6 * np.rad2deg(ras_err) * np.cos(decs)
    decs_err = 3.6e6 * np.rad2deg(decs_err)

    hypothetical = np.linspace(t[0], t[-1], 1000)
    hypothetical_extended = np.linspace(t[0]-1, t[-1]+1, 1000)

    ra_mod, dec_mod = model(hypothetical_extended, ref_t, *result.x)
    ra_mod_long, dec_mod_long = model(hypothetical_extended, ref_t, *result.x)
    ra_mod_wopm, dec_mod_wopm = model_wopm(hypothetical_extended, ref_t, *result.x)
    ra_wopm, dec_wopm = model_subtract_pm(t, ras, decs, *result.x)
    
    # pair observations to points on elipse
    ra_mod_elips, dec_mod_elips = model_wopm(t, ref_t, *result.x)

    params = np.rad2deg(result.x)*3.6e6 
    ra0, dec0, parallax, pm_ra, pm_dec = params
    
    ra0 = np.deg2rad(ra0/3.6e6) 
    dec0 = np.deg2rad(dec0/3.6e6)

    fig = plt.figure(figsize=(15, 10), dpi=120)
    #plt.suptitle(name)

    # Outer grid: 2 blocks (top + bottom)
    outer = gridspec.GridSpec(2, 1, height_ratios=[2, 2], hspace=0.4)

    # Top block (2x2, no spacing)
    top = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[0], hspace=0)

    ax00 = fig.add_subplot(top[0, 0])
    ax10 = fig.add_subplot(top[1, 0], sharex=ax00)

    ax01 = fig.add_subplot(top[0, 1])
    ax11 = fig.add_subplot(top[1, 1], sharex=ax01)

    # Bottom block (normal spacing)
    bottom = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1])

    ax20 = fig.add_subplot(bottom[0, 0])
    ax21 = fig.add_subplot(bottom[0, 1])

    # --- PLOTTING ---

    # With proper motion
    ax00.set_title('Dek och RA med egenrörelse')
    ax00.plot(hypothetical_extended, 3.6e6*np.rad2deg(dec_mod-dec0), color = 'cornflowerblue')
    ax00.errorbar(t, 3.6e6*np.rad2deg(decs-dec0), yerr=decs_err, color='crimson', fmt='o', ms=3, ecolor='black', alpha=0.75)
    ax00.scatter(t, 3.6e6*np.rad2deg(decs-dec0), color='crimson', s=3, zorder = 5)
    ax00.set_ylabel(r'$\Delta\delta$ [mas]')
    ax00.tick_params(labelbottom=False)

    ax10.plot(hypothetical_extended, 3.6e6*np.rad2deg(ra_mod-ra0), color = 'cornflowerblue')
    ax10.errorbar(t, 3.6e6*np.rad2deg(ras-ra0), yerr=ras_err, color='crimson', fmt='o', ms=3, ecolor='black', alpha=0.75)
    ax10.scatter(t, 3.6e6*np.rad2deg(ras-ra0), color='crimson', s=3, zorder = 5)
    ax10.set_ylabel(r'$\Delta\alpha$ [mas]')
    ax10.set_xlabel(r'$\Delta t$ [years]')

    # Without proper motion
    ax01.set_title('Dek och RA utan egenrörelse')
    ax01.plot(hypothetical_extended, 3.6e6*np.rad2deg(dec_mod_wopm-dec0), color = 'cornflowerblue')
    ax01.errorbar(t, 3.6e6*np.rad2deg(dec_wopm-dec0), yerr=decs_err, color='crimson', fmt='o', ms=3, ecolor='black', alpha=0.75)
    ax01.scatter(t, 3.6e6*np.rad2deg(dec_wopm-dec0), color='crimson', s=3, zorder = 5)
    ax01.set_ylabel(r'$\Delta\delta$ [mas]')
    ax01.tick_params(labelbottom=False)

    ax11.plot(hypothetical_extended, 3.6e6*np.rad2deg(ra_mod_wopm-ra0), color = 'cornflowerblue')
    ax11.errorbar(t, 3.6e6*np.rad2deg(ra_wopm-ra0), yerr=ras_err, color='crimson', fmt='o', ms=3, ecolor='black', alpha=0.75)
    ax11.scatter(t, 3.6e6*np.rad2deg(ra_wopm-ra0), color='crimson', s=3, zorder = 5)
    ax11.set_ylabel(r'$\Delta\alpha$ [mas]')
    ax11.set_xlabel(r'$\Delta t$ [years]')

    # Sky plots (separate spacing!)
    ax20.plot(3.6e6*np.rad2deg(ra_mod_long-ra0), 3.6e6*np.rad2deg(dec_mod_long-dec0), color = 'cornflowerblue')
    ax20.errorbar(3.6e6*np.rad2deg(ras-ra0), 3.6e6*np.rad2deg(decs-dec0), xerr=ras_err, yerr=decs_err, color='crimson', fmt='o', ms=3, ecolor='black', alpha=0.4)
    ax20.scatter(3.6e6*np.rad2deg(ras-ra0), 3.6e6*np.rad2deg(decs-dec0), color='crimson', s=3, zorder = 5)
    ax20.axis('equal')
    ax20.invert_xaxis()
    ax20.set_title('Dek mot RA med egenrörelse')
    ax20.set_xlabel(r'$\Delta\alpha$ [mas]')
    ax20.set_ylabel(r'$\Delta\delta$ [mas]')

    ax21.plot(3.6e6*np.rad2deg(ra_mod_wopm-ra0), 3.6e6*np.rad2deg(dec_mod_wopm-dec0), color = 'cornflowerblue', label = 'Modell', linewidth=0.5)
    ax21.errorbar(3.6e6*np.rad2deg(ra_wopm-ra0), 3.6e6*np.rad2deg(dec_wopm-dec0), xerr=ras_err, yerr=decs_err, color='crimson', fmt='o', ms=3, ecolor='black', alpha=0.4)
    ax21.scatter(3.6e6*np.rad2deg(ra_wopm-ra0), 3.6e6*np.rad2deg(dec_wopm-dec0), color='crimson', s=3, label = 'Observationer', zorder = 5)
    if datum:
        # Ta reda på avståden till modellen
        x1=3.6e6*np.rad2deg([dec_wopm-(dec_mod_elips)])
        x2=3.6e6*np.rad2deg([ra_wopm-(ra_mod_elips)])
        dist=np.sqrt(x1**2+x2**2)
        # Ta fram de längsta
        dist0=np.argsort(dist)[:,len(dist[0,:])-1]
        dist1=np.argsort(dist)[:,len(dist[0,:])-2]
        
    for i in range(len(ra_wopm)):
        if datum: 
            if i==dist0 or i==dist1:
                count=0
                # Rätt tid
                tid = aspy.time.Time(ref_t+t[i], format = 'decimalyear')
                tid.format='fits'
                ax21.text(3.6e6*np.rad2deg(ra_wopm[i]-ra0), 3.6e6*np.rad2deg(dec_wopm[i]-dec0), tid)
                # Läs igenom banlistan och se om datumet redan är uppskrivet
                with open(r'banlist_prel.txt', 'r') as fp:
                    lines = fp.readlines()
                    for row in lines:
                        word = f'{tid}'  # Datum som söks efter
                        if row.find(word) != -1:
                            count=count+1
                            print(tid)
                    if count==0:
                        # Om ej uppskrivet så skrivs den in
                        with open("banlist_prel.txt", "a") as f:
                            f.write(f"{tid}\n")
        ax21.plot(3.6e6*np.rad2deg([ra_wopm[i]-ra0, ra_mod_elips[i]-ra0]), 3.6e6*np.rad2deg([dec_wopm[i]-dec0, dec_mod_elips[i]-dec0]), color = 'crimson', alpha=0.25)
    ax21.axis('equal')
    ax21.set_title('Dek mot RA utan egenrörelse')
    ax21.set_xlabel(r'$\Delta\alpha$ [mas]')
    ax21.set_ylabel(r'$\Delta\delta$ [mas]')
    ax21.legend(loc='upper right')
    if save == False:
        plt.show()
    elif save == True:
        plt.savefig(f'Results/starplots/{name}.png')