# Parser
import numpy as np
import scipy as sci
import astropy as aspy
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits


from os import listdir
from os.path import isfile, join
from datetime import datetime


# Returnerar lista av (band (int), tidpunkt (datetime), bild_data (np.ndarray av ra, dec, intensitet))
def parseFits(fits_dir):
    # Hämta lista av fits filer från mapp
    fits_files = [fits_file for fits_file in listdir(fits_dir) if isfile(join(fits_dir, fits_file))]

    observations = list()
    date_format = "%Y-%m-%dT%H:%M:%S.%f"

    # Gå igenom varje fits fil
    for fits_file in fits_files:
        # Öppna fitsfil
        # HDUList (Header Data Unit)
        hdul = fits.open(fits_dir + fits_file)
        #hdul.info() # Printar lite info om fitsfilen

        # Header information
        # Se: http://www.alma.inaf.it/images/ArchiveKeyworkds.pdf
        header = hdul[0].header
        #print(repr(header))

        # Få band
        band = int(fits_file[fits_file.find('B')+1])

        # Få observationsdatum och tid
        date_str = header['DATE-OBS']
        date = datetime.strptime(date_str, date_format)

        # Bilddatan
        intensities = np.squeeze(hdul[0].data)

        width, height = intensities.shape
        x_pixels = np.linspace(0, width, width)
        y_pixels = np.linspace(0, height, height)

        # Right ascension
        ra_ref = header['CRVAL1'] # RA referens i relation till Pixel referensen
        ra_delta = header['CDELT1'] # RA delta för varje pixel
        ra_pixel_ref = header['CRPIX1'] # Pixel positionen där RA referens gäller

        # Declination
        dec_ref = header['CRVAL2']
        dec_delta = header['CDELT2']
        dec_pixel_ref = header['CRPIX2']

        # Konvertera från pixel koordinater till ekvatoriella koordinater
        ras = ra_ref + (x_pixels - ra_pixel_ref) * ra_delta
        decs = dec_ref + (y_pixels - dec_pixel_ref) * dec_delta

        observations.append((band, date, ras, decs, intensities))

    return observations
