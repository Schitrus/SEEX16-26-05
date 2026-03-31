import numpy as np

from datetime import datetime

from utilities import raFormatter, decFormatter



def toAstroDataFile(fname, header, data):
    formatted_data = []
    for data_point in data.T:
        date = data_point[0]
        start = datetime(date.year, 1, 1)
        delta = date - start
        frac_year = date.year + delta.days / (366 if (start.year % 4 == 0) else 365)
        ra = data_point[1]
        dec = data_point[2]
        ra_error = 0.00035 # Dummy value for now
        dec_error = 0.0025 # Dummy value for now
        formatted_data.append((f"{frac_year:.3f}", raFormatter(ra, 0), f"{ra_error:.5f}", decFormatter(dec, 0), f"{dec_error:.5f}"))

    np.savetxt(fname=fname, header="date  RA  RA_error  Dec Dec_error", comments=header + "\n\n# ", X=formatted_data, fmt='%s', delimiter="  ")
