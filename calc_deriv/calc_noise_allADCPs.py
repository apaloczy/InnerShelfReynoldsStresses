# Calculate noise floors for each ADCP from the minimum observed value of the variance.
import numpy as np
from xarray import open_dataset
from ap_tools.utils import m2pydatetime_arr


def get_noise(f, z0, ensminutes=10):
    z0 = dict(z=z0)
    ds = open_dataset(f).sel(z0, method='nearest')
    ds['t'] = m2pydatetime_arr(ds['t'])
    b1, b2, b3, b4, b5 = ds['b1'], ds['b2'], ds['b3'], ds['b4'], ds['b5']
    b1var = (b1 - b1.mean('t'))**2
    b2var = (b2 - b2.mean('t'))**2
    b3var = (b3 - b3.mean('t'))**2
    b4var = (b4 - b4.mean('t'))**2
    b5var = (b5 - b5.mean('t'))**2
    b1var = b1var.resample(dict(t="%d min"%ensminutes)).mean('t')
    b2var = b2var.resample(dict(t="%d min"%ensminutes)).mean('t')
    b3var = b3var.resample(dict(t="%d min"%ensminutes)).mean('t')
    b4var = b4var.resample(dict(t="%d min"%ensminutes)).mean('t')
    b5var = b5var.resample(dict(t="%d min"%ensminutes)).mean('t')
    Nb1 = np.sqrt(b1var.min().values.flatten()[0])
    Nb2 = np.sqrt(b2var.min().values.flatten()[0])
    Nb3 = np.sqrt(b3var.min().values.flatten()[0])
    Nb4 = np.sqrt(b4var.min().values.flatten()[0])
    Nb5 = np.sqrt(b5var.min().values.flatten()[0])
    Nb = 0.25*(Nb1 + Nb2 + Nb3 + Nb4)

    return Nb, Nb5 # [m/s], an rms noise.


#---
head = '/home/andre/phd/data/innershelfdri_moorings'
z025 = 8
z040 = 15

# OC25SA.
f1 = head + '/OC25SA-A/deployment1/OC25SA-Ad1beam.nc'
f2 = head + '/OC25SA-A/deployment2/OC25SA-Ad2beam.nc'
Nb_depl1, Nb5_depl1 = get_noise(f1, z025)
Nb_depl2, Nb5_depl2 = get_noise(f2, z025)
Nb = (Nb_depl1 + Nb_depl2)/2
Nb5 = (Nb5_depl1 + Nb5_depl2)/2
np.savez("noise_OC25SA.npz", Nb=Nb, Nb5=Nb5)
print("")
print("OC25SA: Nb = %.4f, Nb5 = %.4f m/s"%(Nb, Nb5))

# OC25SB.
f2 = head + '/OC25SB-A/deployment2/OC25SB-Ad2beam.nc'
Nb, Nb5 = get_noise(f2, z025)
np.savez("noise_OC25SB.npz", Nb=Nb, Nb5=Nb5)
print("")
print("OC25SB: Nb = %.4f, Nb5 = %.4f m/s"%(Nb, Nb5))

# OC40N.
f1 = head + '/OC40N-A/deployment1/OC40N-Ad1beam.nc'
Nb, Nb5 = get_noise(f1, z040)
np.savez("noise_OC40N.npz", Nb=Nb, Nb5=Nb5)
print("")
print("OC40N: Nb = %.4f, Nb5 = %.4f m/s"%(Nb, Nb5))

# OC40S.
f1 = head + '/OC40S-A/deployment1/OC40S-Ad1beam.nc'
f2 = head + '/OC40S-A/deployment2/OC40S-Ad2beam.nc'
Nb_depl1, Nb5_depl1 = get_noise(f1, z040)
Nb_depl2, Nb5_depl2 = get_noise(f2, z040)
Nb = (Nb_depl1 + Nb_depl2)/2
Nb5 = (Nb5_depl1 + Nb5_depl2)/2
np.savez("noise_OC40S.npz", Nb=Nb, Nb5=Nb5)
print("")
print("OC40S: Nb = %.4f, Nb5 = %.4f m/s"%(Nb, Nb5))
