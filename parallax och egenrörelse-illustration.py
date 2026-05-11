import numpy as np
import matplotlib.pyplot as plt

import parallax_fitting as pf

t = np.linspace(0, 10, 1000)
ref_t = 2000
alpha_0 = np.deg2rad(69.19032363779596) # mas to radians
delta_0 = np.deg2rad(-62.076981747826686) # mas to radians
parallax = np.deg2rad(18.31/3.6e6) # mas to radians
pm_ra = np.deg2rad(-69.36/3.6e6) # mas/yr to radians
pm_dec = np.deg2rad(-75.78/3.6e6) # mas/yr to radians
params = [alpha_0, delta_0, parallax, pm_ra, pm_dec]
result = pf.model(t, ref_t, *params)
result_wopm = pf.model_wopm(t, ref_t, *params)

fig, ax = plt.subplots(1,2, figsize=(12, 6), dpi=120)
ax[0].set_title('Med egenrörelse')
ax[0].plot(3.6e6*np.rad2deg(result[0] - alpha_0), 3.6e6*np.rad2deg(result[1] - delta_0), color = 'cornflowerblue', label=rf'$\mu_\alpha$ = {np.rad2deg(pm_ra)*3.6e6:.2f} mas/yr, $\mu_\delta$ = {np.rad2deg(pm_dec)*3.6e6:.2f} mas/yr'.replace(".", ","))
ax[0].invert_xaxis()
ax[0].set_xlabel(rf'$\Delta\alpha$ (mas)')
ax[0].set_ylabel(rf'$\Delta\delta$ (mas)')
ax[0].legend(loc = 'upper right')
ax[1].set_title('Utan egenrörelse')
ax[1].plot(3.6e6*np.rad2deg(result_wopm[0] - alpha_0), 3.6e6*np.rad2deg(result_wopm[1] - delta_0), color = 'cornflowerblue', label = rf'$\varpi$ = {np.rad2deg(parallax)*3.6e6:.2f} mas'.replace(".", ","))
ax[1].set_xlabel(rf'$\Delta\alpha$ (mas)')
ax[1].set_ylabel(rf'$\Delta\delta$ (mas)')
ax[1].legend(loc = 'upper right')
plt.savefig('Exempelgrafer-egenrörelse-parallax.png')  
#plt.show()

