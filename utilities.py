# Formatterar deklinationen till string format
def decFormatter(dec, pos): 
    dec_deg = int(dec)
    dec_arcmin = int(60 * (dec - dec_deg))
    dec_arcsec = 60 * (60 * (dec - dec_deg) - dec_arcmin)

    #return f"{dec_deg}°{abs(dec_arcmin)}'{abs(dec_arcsec):.3f}''"
    return f"{dec_deg}:{abs(dec_arcmin)}:{abs(dec_arcsec):.6f}"

def decDeformatter(dec):
    angle, arcmin, arcsec = dec.split(':')
    return float(angle) + float(arcmin)/60 + float(arcsec)/3600

# Formatterar rektascension till string format
def raFormatter(ra, pos): 
    ra /= 15
    ra_hour = int(ra)
    ra_min = int(60 * (ra - ra_hour))
    ra_sec = 60 * (60 * (ra - ra_hour) - ra_min)

    #return f"{ra_hour}h{abs(ra_min)}m{abs(ra_sec):.3f}s"
    return f"{ra_hour:02d}:{abs(ra_min)}:{abs(ra_sec):.6f}"

def raDeformatter(ra):
    hour, arcmin, arcsec = ra.split(':')
    return 15*(float(hour) + float(arcmin)/60 + float(arcsec)/3600)