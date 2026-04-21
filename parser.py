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
from PIL import Image

#läs banlist
with open(r'banlist.txt', 'r') as fp:
    lines = fp.readlines()

# Returnerar lista av (band (int), tidpunkt (datetime), bild_data (np.ndarray av ra, dec, intensitet))
def parseFits(fits_dir, user_file = None, banlist=False):
    # Hämta lista av fits filer från mapp
    fits_files = [fits_file for fits_file in listdir(fits_dir) if isfile(join(fits_dir, fits_file))]
    observations = list()
    date_format = "%Y-%m-%dT%H:%M:%S.%f"

    if user_file:
        user_data = np.loadtxt(fname='user_data/' + user_file, skiprows= 3, unpack=True, encoding='latin1', 
                                                    dtype=[("time", '<U20'), ("ra_pixels",float), ("dec_pixels", float)])
        user_name, user_ra_pixels, user_dec_pixels = user_data

    # Gå igenom varje fits fil
    for fits_file in fits_files:
        next=False
        # Kolla om nuvarande fil är i banlist
        for i in range(len(lines)):
            if f"{fits_file}\n"==f"{lines[i]}":
                next=True
        # Om den är med så gå vidare till nästa
        if banlist and next:
            continue
        # Annars fortsätt som vanligt
        else:
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
            im = Image.fromarray((255*(intensities - np.min(intensities))/(np.max(intensities)-np.min(intensities))).astype(dtype=np.uint8))
            im.save(f'star_images/{fits_file}.png')

        

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

            if user_file:
                idx = np.where(user_name == fits_file.split('.')[0])
                if idx[0].size == 0:
                    continue
                u_ra_pixel, u_dec_pixel = user_ra_pixels[idx][0], user_dec_pixels[idx][0]

                user_ra = ra_ref + (u_ra_pixel - ra_pixel_ref) * ra_delta
                user_dec = dec_ref + (u_dec_pixel - dec_pixel_ref) * dec_delta

                observations.append((band, date, ras, decs, intensities, np.abs(ra_delta), np.abs(dec_delta), user_ra, user_dec))
            else:
                observations.append((band, date, ras, decs, intensities))

    return observations
