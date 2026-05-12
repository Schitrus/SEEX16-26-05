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
        ra_error = aspy.coordinates.Angle(data_point[3]*aspy.units.degree).to_string(precision=16)
        dec_error = aspy.coordinates.Angle(data_point[4]*aspy.units.degree).to_string(precision=16)
        formatted_data.append((f"{decimal_year:.3f}", ra, ra_error, dec, dec_error))

    np.savetxt(fname=fname, header="date  RA  RA_error  Dec Dec_error", comments=header + "\n\n# ", X=formatted_data, fmt='%s', delimiter="  ")

def fromAstroData(fname):
    with open(fname, "r") as f:
        first_line = f.readline().strip()

    name = first_line.split("=", 1)[1].strip()
    formatter = {0: lambda ts: aspy.time.Time(float(ts), format='decimalyear'), 
                 1: lambda ras: aspy.coordinates.Angle(ras, unit=aspy.units.deg).degree,
                 2: lambda ra_errors: aspy.coordinates.Angle(ra_errors, unit=aspy.units.deg).degree,
                 3: lambda decs: aspy.coordinates.Angle(decs, unit=aspy.units.deg).degree,
                 4: lambda dec_errors: aspy.coordinates.Angle(dec_errors, unit=aspy.units.deg).degree}
    data = np.loadtxt(fname=fname, skiprows= 3, unpack=True, converters=formatter,encoding='latin1', 
                                                     dtype=[("time", aspy.time.Time), ("ras",float), ("ra_errs", float), ("decs", float), ("dec_errs", float)])
    ts, ras, ra_errors, decs, dec_errors = data
    return name, ts, ras, decs, ra_errors, dec_errors