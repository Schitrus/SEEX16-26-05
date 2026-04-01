import numpy as np

from datetime import datetime

import astropy as aspy

from utilities import raFormatter, decFormatter, raDeformatter, decDeformatter

def toAstroData(fname, header, data):
    formatted_data = []
    for data_point in data.T:
        date = aspy.time.Time(data_point[0])
        decimal_year = date.decimalyear
        
        ra = aspy.coordinates.Angle(data_point[1]*aspy.units.degree).to_string(precision=16)
        dec = aspy.coordinates.Angle(data_point[2]*aspy.units.degree).to_string(precision=16)
        ra_error = 0.00035 # Dummy value for now
        dec_error = 0.0025 # Dummy value for now
        formatted_data.append((f"{decimal_year:.3f}", ra, f"{ra_error:.5f}", dec, f"{dec_error:.5f}"))

    np.savetxt(fname=fname, header="date  RA  RA_error  Dec Dec_error", comments=header + "\n\n# ", X=formatted_data, fmt='%s', delimiter="  ")

def fromAstroData(fname):
    formatter = {0 : lambda ts: aspy.time.Time(float(ts), format='decimalyear'), 
                 1: lambda ras: aspy.coordinates.Angle(ras).degree, 
                 2: lambda ra_errors: 15*float(ra_errors),
                 3: lambda decs: aspy.coordinates.Angle(decs).degree,
                 4: lambda dec_errors: float(dec_errors)}
    data = np.loadtxt(fname=fname, skiprows= 3, unpack=True, converters=formatter, 
                                                     dtype=[("time", aspy.time.Time), ("ras",float), ("ra_errs", float), ("decs", float), ("dec_errs", float)])
    ts, ras, ra_errors, decs, dec_errors = data
    return ts, ras, decs