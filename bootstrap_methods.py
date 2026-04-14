import numpy as np

from parallax_fitting import model, fit_model, residuals, bary_coords



# Function for using bootstrap to get uncertainties in data
def bootstrap(t, ref_t, ras, decs, result, bounds, ras_err, decs_err, n_boot=1000):
    
    ra_model, dec_model = model(t, ref_t, *result.x)

    res_ra = ras - ra_model
    res_dec = decs - dec_model

    n = len(t)

    param_samples = []
    
    for i in range(n_boot):
        idx = np.random.randint(0,n,n)

        ra_boot = ra_model + res_ra[idx]
        dec_boot = dec_model + res_dec[idx]

        ras_err_boot = ras_err[idx]
        decs_err_boot = decs_err[idx]
        try:
            res_boot = fit_model(t, ref_t, ra_boot, dec_boot, result.x, bounds, ras_err_boot, decs_err_boot)

            param_samples.append(res_boot.x)
        except Exception:
            continue

    print(f'Bootstrap success: {len(param_samples)}/{n_boot}')
    
    return np.array(param_samples)

# Bootstrap with leave-one-out method
def leave_one_out(t, ref_t, ras, decs, result, bounds, ras_err, decs_err):
    
    n = len(t)
    param_samples = []

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False

        try:
            res = fit_model(
                t[mask], ref_t,
                ras[mask], decs[mask],
                result.x, bounds,
                ras_err[mask], decs_err[mask]
            )

            param_samples.append(res.x)

        except Exception:
            continue

    print(f'LOO success: {len(param_samples)}/{n}')
    
    return np.array(param_samples)

# Bootstrap with leave-one-epoch-out method
def leave_one_epoch_out(t, ref_t, ras, decs, result, bounds, ras_err, decs_err):
    
    unique_times = np.unique(t)
    param_samples = []

    for time in unique_times:
        mask = t != time

        try:
            res = fit_model(
                t[mask], ref_t,
                ras[mask], decs[mask],
                result.x, bounds,
                ras_err[mask], decs_err[mask]
            )

            param_samples.append(res.x)

        except Exception:
            continue

    print(f'LOEO success: {len(param_samples)}/{len(unique_times)}')
    
    return np.array(param_samples)