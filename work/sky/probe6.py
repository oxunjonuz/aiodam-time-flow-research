import numpy as np, math
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u
h=fits.open('LALInference_skymap.fits')
prob=np.asarray(h[1].data['PROB']).reshape(-1).astype(float)
nside=512
hp=HEALPix(nside=nside, order='nested', frame='icrs')
pix=np.arange(len(prob))
lon,lat=hp.healpix_to_lonlat(pix)
ra=lon.to_value(u.rad); dec=lat.to_value(u.rad)
n=np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)],axis=1)
C=299792458.0
VERT={"H1":np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6]),
      "L1":np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])}
delay=(VERT["H1"]-VERT["L1"])@n.T/C*1e3  # ms, positive = H1 later
w=prob/prob.sum()
print("delay over full sky: min %.3f max %.3f"%(delay.min(),delay.max()))
print("prob-weighted mean delay %.4f ms"%float((delay*w).sum()))
order=np.argsort(prob)[::-1]
cum=np.cumsum(prob[order])
n90=int(np.searchsorted(cum,0.9))+1
sel=order[:n90]
print("90%% region: delay min %.3f max %.3f mean %.4f"%(delay[sel].min(),delay[sel].max(),float((delay[sel]*prob[sel]).sum()/prob[sel].sum())))
print("top-20 pixels delay:", np.round(delay[order[:20]],3))
print("top-20 RA/Dec:", np.round(np.degrees(ra[order[:20]]),2), np.round(np.degrees(dec[order[:20]]),2))
# how much probability sits within +-0.5 ms of 6.9?
m=np.abs(delay-6.9)<0.5
print("prob within 0.5ms of 6.9ms: %.4f"%prob[m].sum())
m2=np.abs(delay-1.25)<0.5
print("prob within 0.5ms of 1.25ms: %.4f"%prob[m2].sum())
