import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from xarray import DataArray, Dataset, open_dataset

from dewaveADCP.stress import rot_uwvw
from dewaveADCP.AdaptiveFiltering import bvelAF
from dewaveADCP.beam2earth import getmskfish, binmap, binmap5

import warnings
warnings.filterwarnings("ignore")

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)

def m2pydatetime_arr(mdatenum_arr):
    mdatenum_arr = np.array(mdatenum_arr)
    timedt = timedelta(days=366)
    pydt = []
    for mdt in mdatenum_arr.tolist():
        d = datetime.fromordinal(int(mdt))
        dfrac = timedelta(days=mdt%1) - timedt
        pydt.append(d + dfrac)

    return np.array(pydt)


def sind(x):
    return np.sin(x*np.pi/180)


def cosd(x):
    return np.cos(x*np.pi/180)


def bvar(b):
    """
    Calculates the time series of profiles of variance of along-beam velocity
    from the along-beam velocities 'b'.
    """
    return (b - np.nanmean(b, axis=1)[:, np.newaxis])**2


def bvarm(b):
    """
    Calculates the time-averaged profiles of variance of along-beam velocity
    from the along-beam velocities 'b'.
    """
    return np.nanmean((b - np.nanmean(b, axis=1)[:, np.newaxis])**2, axis=1)


def uwrs4(b1var, b2var, theta):
    Sth, Cth = sind(theta), cosd(theta)
    b2mb1 = b2var - b1var
    coeff = 1/(4*Sth*Cth)

    # Dewey & Stringer (2007)'s equation (32).
    uw = -coeff*b2mb1

    return uw


def vwrs4(b3var, b4var, theta):
    Sth, Cth = sind(theta), cosd(theta)
    b4mb3 = b4var - b3var
    coeff = 1/(4*Sth*Cth)

    # Dewey & Stringer (2007)'s equation (33).
    vw = -coeff*b4mb3

    return vw


def uwrs5(b1, b2, b5, uv, theta, phi2, phi3, variances=True, Nb=0, Nb5=0):
    if variances:
        b1var, b2var, b5var = b1, b2, b5
    else:
        b1var, b2var, b5var = map(bvar, (b1, b2, b5))

    Sth, Cth = sind(theta), cosd(theta)
    S2 = Sth**2
    S6C2 = (Sth**6)*(Cth**2)
    S5C1 = (Sth**5)*(Cth)
    S4C2 = (Sth**4)*(Cth**2)

    phi2, phi3 = phi2*d2r, phi3*d2r
    b2mb1 = b2var - b1var
    b2pb1 = b2var + b1var
    coeff = -1/(4*S6C2)

    # Dewey & Stringer (2007)'s equation (132).
    uw = coeff*(S5C1*b2mb1 + 2*S4C2*phi3*b2pb1 - 4*S4C2*phi3*b5var - 4*S6C2*phi2*uv)

    phi3ens = np.nanmean(phi3)
    bias = phi3ens*(Nb5**2 - Nb**2)/S2

    return np.nanmean(uw, axis=1) - bias


def vwrs5(b3, b4, b5, uv, theta, phi2, phi3, variances=True, Nb=0, Nb5=0):
    if variances:
        b3var, b4var, b5var = b3, b4, b5
    else:
        b3var, b4var, b5var = map(bvar, (b3, b4, b5))

    Sth, Cth = sind(theta), cosd(theta)
    S2, C2 = Sth**2, Cth**2
    S6C2 = (Sth**6)*(Cth**2)
    S5C1 = (Sth**5)*(Cth)
    S4C2 = (Sth**4)*(Cth**2)
    S4C4 = (Sth*Cth)**4

    phi2, phi3 = phi2*d2r, phi3*d2r
    b4mb3 = b4var - b3var
    b4pb3 = b4var + b3var
    coeff = -1/(4*S6C2)

    # Dewey & Stringer (2007)'s equation (133).
    vw = coeff*(S5C1*b4mb3 - 2*S4C2*phi3*b4pb3 + 4*S4C4*phi3*b5var + 4*S6C2*phi2*b5var + 4*S6C2*phi3*uv)

    phi2ens, phi3ens = np.nanmean(phi2), np.nanmean(phi3)
    bias = phi3ens*(Nb**2 - C2*Nb5**2)/S2 - phi2ens*Nb5**2

    return np.nanmean(vw, axis=1) - bias



def tke5(b1, b2, b3, b4, b5, theta, phi3, variances=True, Nb=0, Nb5=0):
    if variances:
        b1var, b2var, b3var, b4var, b5var = b1, b2, b3, b4, b5
    else:
        b1var, b2var, b3var, b4var, b5var = map(bvar, (b1, b2, b3, b4, b5))

    Sth, Cth = sind(theta), cosd(theta)
    S2 = Sth**2
    C2 = Cth**2
    cotth = Cth/Sth

    phi3 = phi3*d2r
    b2mb1 = b2var - b1var
    b1234 = b1var + b2var + b3var + b4var
    coeff = 1/(4*S2)

    # D&S Equation 134.
    q2 = coeff*(b1234 - 2*(2*C2 - S2)*b5var - (cotth - 1)*phi3*b2mb1) # q^2/2, not q^2.

    bias = coeff*(4*Nb**2 - 2*(2*C2 - S2)*Nb5**2)

    return np.nanmean(q2, axis=1) - bias


def aniso_ratio(b1, b2, b3, b4, b5, theta, phi2, phi3, variances=True):
    if variances:
        b1var, b2var, b3var, b4var, b5var = b1, b2, b3, b4, b5
    else:
        b1var, b2var, b3var, b4var, b5var = map(bvar, (b1, b2, b3, b4, b5))

    Sth, Cth = sind(theta), cosd(theta)
    S2 = Sth**2
    C2 = Cth**2
    csc2th = 1/sind(2*theta)
    Tth = Sth/Cth

    phi2, phi3 = phi2*d2r, phi3*d2r
    b2mb1 = b2var - b1var
    b4mb3 = b4var - b3var
    b1234 = b1var + b2var + b3var + b4var
    Fth = Tth*phi2*b4mb3 + (1 - 2*csc2th)*phi3*b2mb1

    # D&S Equation 135.
    num = 2*S2*b5var + Tth*phi3*b2mb1 - Tth*phi2*b4mb3
    den = b1234 - 4*C2*b5var + Fth

    alpha = num/den

    return np.nanmean(alpha, axis=1)


#---
plt.close('all')

###
dtens_min = 10
###
BINMAP = None # 'linear' or 'nearest' or None.
freq = 1 # [Hz]
max_badfrac = 0.5
MSKFISH = False
MSKFISH_THRESH = 10
CORRECT_BIAS = True
###

Lw = 8 # Used in the Adaptive Filtering Method.
theta = 25
N = int(freq*60*dtens_min)

# Load dataset.
# fname = '/home/andre/phd/data/innershelfdri_moorings/OC40S-A/deployment1/OC40S-Ad1beam.nc'; fname_enu = '/home/andre/phd/data/innershelfdri_moorings/OC40S-A/deployment1/OC40S-Ad1enu.nc'; sep = 2

fname = '/home/andre/phd/data/innershelfdri_moorings/OC40S-A/deployment2/OC40S-Ad2beam.nc'; fname_enu = '/home/andre/phd/data/innershelfdri_moorings/OC40S-A/deployment2/OC40S-Ad2enu.nc'; sep = 2

# fname = '/home/andre/phd/data/innershelfdri_moorings/OC40N-A/deployment1/OC40N-Ad1beam.nc'; fname_enu = '/home/andre/phd/data/innershelfdri_moorings/OC40N-A/deployment1/OC40N-Ad1enu.nc'; sep = 2

# fname = '/home/andre/phd/data/innershelfdri_moorings/OC25SA-A/deployment1/OC25SA-Ad1beam.nc'; fname_enu = '/home/andre/phd/data/innershelfdri_moorings/OC25SA-A/deployment1/OC25SA-Ad1enu.nc'; sep = 4

# fname = '/home/andre/phd/data/innershelfdri_moorings/OC25SA-A/deployment2/OC25SA-Ad2beam.nc'; fname_enu = '/home/andre/phd/data/innershelfdri_moorings/OC25SA-A/deployment2/OC25SA-Ad2enu.nc'; sep = 4

# fname = '/home/andre/phd/data/innershelfdri_moorings/OC25SB-A/deployment2/OC25SB-Ad2beam.nc'; fname_enu = '/home/andre/phd/data/innershelfdri_moorings/OC25SB-A/deployment2/OC25SB-Ad2enu.nc'; sep = 4

# Get noise floors of Janus and vertical beam for each ADCP.
if CORRECT_BIAS:
    moor = fname.split('/')[-1].split('-')[0]
    dnoise = np.load("noise_%s.npz"%moor)
    Nb, Nb5 = dnoise["Nb"].flatten()[0], dnoise["Nb5"].flatten()[0]
else:
    Nb, Nb5 = 0, 0

tl, tr = '2017-01-12 00:00:00', '2017-12-12 03:00:00'
# tl, tr = '2017-09-10 00:00:00', '2017-09-12 03:00:00'
# tl, tr = '2017-10-15 00:00:00', '2017-10-17 00:00:00'
# tl, tr = '2017-10-22 01:00:00', '2017-10-22 06:10:00'

kwbvelAF = dict(sep=sep, Lw=Lw, verbose=False, max_badfrac=max_badfrac)
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

fname_rs_out = fname.split('/')[-1].split('.')[0].replace('beam', 'rs5')
fname_rs_out = fname_rs_out + '-%dmin'%dtens_min + '_ensembles.nc'
if BINMAP is None:
    fname_rs_out = fname_rs_out.split('.')[0] + '-nobinmap.nc'
else:
    fname_rs_out = fname_rs_out.split('.')[0] + '-%sbinmap.nc'%BINMAP

if CORRECT_BIAS:
    fname_rs_out = fname_rs_out.split('.')[0] + '-bias-corrected.nc'
else:
    fname_rs_out = fname_rs_out.split('.')[0] + '-not-bias-corrected.nc'

if True:
    fname_rs_out = "_sample.".join(fname_rs_out.split("."))

ds = open_dataset(fname)
ds['t'] = m2pydatetime_arr(ds['t'].values)

ds = ds.sel(dict(t=slice(tl, tr)))

tj, rj, zj, b1, b2, b3, b4, b5, hdng, ptch, roll, p = ds['t'].values, ds['r'].values, ds['z'].values, ds['b1'].T, ds['b2'].T, ds['b3'].T, ds['b4'].T, ds['b5'].T, ds['hdng'], ds['ptch'], ds['roll'], ds['p']
amp1, amp2, amp3, amp4 = ds['amp1'].T, ds['amp2'].T, ds['amp3'].T, ds['amp4'].T

dsenu = open_dataset(fname_enu)
dsenu['t'] = m2pydatetime_arr(dsenu['t'].values)
dsenu = dsenu.sel(dict(t=slice(tl, tr)))

u, v, w5 = dsenu['u'].T, dsenu['v'].T, dsenu['w5'].T
h = zj[-1]

nz, nt = b1.shape
nens = int(nt//N) # Number of ensembles in stress time series.

uw = np.empty((nz, nens))*np.nan
vw = uw.copy()
uwraw, vwraw = vw.copy(), vw.copy()
uwrawnotilt, vwrawnotilt = vw.copy(), vw.copy()
uwnotilt, vwnotilt = vw.copy(), vw.copy()
uw, vw = vw.copy(), vw.copy()
tke, aniso = uw.copy(), uw.copy()
uens, vens, upvpens = vw.copy(), vw.copy(), vw.copy()
w5ens = uw.copy()
urms, vrms, w5rms = vw.copy(), vw.copy(), vw.copy()
pens = np.empty(nens)*np.nan
hdngens, ptchens, rollens = pens.copy(), pens.copy(), pens.copy()

bad = np.empty(nz)*np.nan
tjens = np.array([], dtype=np.datetime64)

d2r = np.pi/180
Sth = np.sin(theta*d2r)
n = 0
il = 0
ir = il + N
while ir<nt:
    print("Ensemble ",n+1," / ",nens)
    # ===Adaptive Filtering Method===
    hdngo = hdng[il:ir].values
    ptcho = ptch[il:ir].values
    rollo = roll[il:ir].values
    po = p[il:ir].values

    tjo = tj[il:ir]
    b1o = b1[:, il:ir].values
    b2o = b2[:, il:ir].values
    b3o = b3[:, il:ir].values
    b4o = b4[:, il:ir].values
    b5o = b5[:, il:ir].values
    w5o = w5[:, il:ir]

    if BINMAP is not None:
        thetar = np.radians(theta)
        ptchor = np.radians(ptch)
        rollor = np.radians(roll)
        b1o, b2o, b3o, b4o, b5o = binmap5(b1o, b2o, b3o, b4o, b5o, rj, thetar, ptchor, rollor, how=BINMAP)

    if MSKFISH:
        if BINMAP is not None:
            amp1, amp2, amp3, amp4 = binmap(amp1, amp2, amp3, amp4, rj, thetar, ptchor, rollor, how=BINMAP)
        amp1o = np.nanmean(amp1[:, il:ir], axis=1, keepdims=True)
        amp2o = np.nanmean(amp2[:, il:ir], axis=1, keepdims=True)
        amp3o = np.nanmean(amp3[:, il:ir], axis=1, keepdims=True)
        amp4o = np.nanmean(amp4[:, il:ir], axis=1, keepdims=True)

        # Calculate mask to eliminate unrealistic vertical stripes.
        msk1, msk2, msk3, msk4 = getmskfish(amp1o, amp2o, amp3o, amp4o, threshold=MSKFISH_THRESH)
        msk = msk1*msk2*msk3*msk4
        msk = msk.flatten()
    else:
        msk = 1.0

    bb1 = bvelAF(b1o, **kwbvelAF)
    bb2 = bvelAF(b2o, **kwbvelAF)
    bb3 = bvelAF(b3o, **kwbvelAF)
    bb4 = bvelAF(b4o, **kwbvelAF)
    bb5 = bvelAF(b5o, **kwbvelAF)

    uwnotilt[:, n] = uwrs4(bvarm(bb1), bvarm(bb2), theta)*msk
    vwnotilt[:, n] = vwrs4(bvarm(bb3), bvarm(bb4), theta)*msk

    ua, va = (-bb1 + bb2)*Sth/2, (-bb3 + bb4)*Sth/2
    uv = (ua - np.nanmean(ua, axis=1, keepdims=True))*(va - np.nanmean(va, axis=1, keepdims=True))
    uw[:, n] = uwrs5(bb1, bb2, bb5, uv, theta, ptcho, rollo, variances=False, Nb=Nb, Nb5=Nb5)*msk
    vw[:, n] = vwrs5(bb3, bb4, bb5, uv, theta, ptcho, rollo, variances=False, Nb=Nb, Nb5=Nb5)*msk

    tke[:, n] = tke5(bb1, bb2, bb3, bb4, bb5, theta, rollo, variances=False, Nb=Nb, Nb5=Nb5)*msk
    # aniso[:, n] = aniso_ratio(bb1, bb2, bb3, bb4, bb5, theta, ptcho, rollo, variances=False)*msk

    # Non-corrected RS.
    uwrawnotilt[:, n] = uwrs4(bvarm(b1o), bvarm(b2o), theta)*msk
    vwrawnotilt[:, n] = vwrs4(bvarm(b3o), bvarm(b4o), theta)*msk

    ua, va = (-b1o + b2o)*Sth/2, (-b3o + b4o)*Sth/2
    uvraw = (ua - np.nanmean(ua, axis=1, keepdims=True))*(va - np.nanmean(va, axis=1, keepdims=True))
    uwraw[:, n] = uwrs5(b1o, b2o, b5o, uvraw, theta, ptcho, rollo, variances=False, Nb=Nb, Nb5=Nb5)*msk
    vwraw[:, n] = vwrs5(b3o, b4o, b5o, uvraw, theta, ptcho, rollo, variances=False, Nb=Nb, Nb5=Nb5)*msk

    # Ensemble-averaged ENU velocities and other quantities.
    uo, vo = u[:, il:ir], v[:, il:ir]
    uens[:, n] = np.nanmean(uo, axis=1) # EARTH coordinates.
    vens[:, n] = np.nanmean(vo, axis=1) # EARTH coordinates.
    w5ens[:, n] = np.nanmean(w5o, axis=1)

    urms[:, n] = np.nanstd(uo, axis=1)
    vrms[:, n] = np.nanstd(vo, axis=1)
    w5rms[:, n] = np.nanstd(w5o, axis=1)

    hdngens[n] = np.nanmean(hdngo)
    ptchens[n] = np.nanmean(ptcho)
    rollens[n] = np.nanmean(rollo)
    pens[n] = np.nanmean(po)

    tjens = np.append(tjens, tj[il+N//2])

    il = ir
    ir += N
    if ir>nt: # Last ensemble might have less pings than the others.
        ir = nt
        nscrap = (nt - il)/60
        print("Last ensemble was %.2f min long."%nscrap)
    n += 1

uwnotilt, vwnotilt = rot_uwvw(uwnotilt, vwnotilt, hdngens)
uw, vw = rot_uwvw(uw, vw, hdngens)
uwrawnotilt, vwrawnotilt = rot_uwvw(uwrawnotilt, vwrawnotilt, hdngens)
uwraw, vwraw = rot_uwvw(uwraw, vwraw, hdngens)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Rotate to across-along-isobath coordinate system.
# Rotation angle at 25 m isobath: -16.72 deg
# Rotation angle at 40 m isobath: -5.90 deg
if "OC25" in fname:
    ang_isob = -16.72
elif "OC40" in fname:
    ang_isob = -5.90

uwnotilt, vwnotilt = rot_uwvw(uwnotilt, vwnotilt, ang_isob)
uw, vw = rot_uwvw(uw, vw, ang_isob)
uwrawnotilt, vwrawnotilt = rot_uwvw(uwrawnotilt, vwrawnotilt, ang_isob)
uwraw, vwraw = rot_uwvw(uwraw, vwraw, ang_isob)
uens, vens = rot_uwvw(uens, vens, ang_isob)
urms, vrms = rot_uwvw(urms, vrms, ang_isob)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Calculate shear.
dz = np.gradient(zj)[:,np.newaxis]
uzens = np.gradient(uens, axis=0)/dz # EARTH coordinates.
vzens = np.gradient(vens, axis=0)/dz # EARTH coordinates.

#### Save stresses and uncertainties in tilt correction terms.

dims = ('z', 't')
coords = dict(t=tjens, z=zj, r=('z', rj))
dims1 = 't'
coords1 = dict(t=tjens)

uwrawnotilt = DataArray(uwrawnotilt, coords=coords, dims=dims)
vwrawnotilt = DataArray(vwrawnotilt, coords=coords, dims=dims)
uwraw = DataArray(uwraw, coords=coords, dims=dims)
vwraw = DataArray(vwraw, coords=coords, dims=dims)
uwnotilt = DataArray(uwnotilt, coords=coords, dims=dims)
vwnotilt = DataArray(vwnotilt, coords=coords, dims=dims)
uw = DataArray(uw, coords=coords, dims=dims)
vw = DataArray(vw, coords=coords, dims=dims)
tke = DataArray(tke, coords=coords, dims=dims)
# aniso = DataArray(aniso, coords=coords, dims=dims)
u = DataArray(uens, coords=coords, dims=dims)
v = DataArray(vens, coords=coords, dims=dims)
w5 = DataArray(w5ens, coords=coords, dims=dims)
urms = DataArray(urms, coords=coords, dims=dims)
vrms = DataArray(vrms, coords=coords, dims=dims)
w5rms = DataArray(w5rms, coords=coords, dims=dims)
uz = DataArray(uzens, coords=coords, dims=dims)
vz = DataArray(vzens, coords=coords, dims=dims)

hdng = DataArray(hdngens, coords=coords1, dims=dims1)
ptch = DataArray(ptchens, coords=coords1, dims=dims1)
roll = DataArray(rollens, coords=coords1, dims=dims1)
p = DataArray(pens, coords=coords1, dims=dims1)

# data_vars = dict(uwrawnotilt=uwrawnotilt, vwrawnotilt=vwrawnotilt, uwraw=uwraw, vwraw=vwraw, uwnotilt=uwnotilt, vwnotilt=vwnotilt, uw=uw, vw=vw, tke=tke, aniso=aniso, u=u, v=v, w5=w5, urms=urms, vrms=vrms, w5rms=w5rms, uz=uz, vz=vz, hdng=hdng, ptch=ptch, roll=roll, p=p)

data_vars = dict(uwrawnotilt=uwrawnotilt, vwrawnotilt=vwrawnotilt, uwraw=uwraw, vwraw=vwraw, uwnotilt=uwnotilt, vwnotilt=vwnotilt, uw=uw, vw=vw, tke=tke, u=u, v=v, w5=w5, urms=urms, vrms=vrms, w5rms=w5rms, uz=uz, vz=vz, hdng=hdng, ptch=ptch, roll=roll, p=p)

Dataset(data_vars=data_vars, coords=coords).to_netcdf(fname_rs_out)
