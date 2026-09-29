import numpy as np, math
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u
h=fits.open('LALInference_skymap.fits')
raw=np.asarray(h[1].data['PROB']).astype(float)
print("raw shape",raw.shape)
C=299792458.0
H1=np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6])
L1=np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])
hp_n=HEALPix(nside=512, order='nested', frame='icrs')
hp_r=HEALPix(nside=512, order='ring', frame='icrs')
pix=np.arange(3145728)
def delay_for(hp):
    lon,lat=hp.healpix_to_lonlat(pix)
    ra=lon.to_value(u.rad); dec=lat.to_value(u.rad)
    n=np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)],axis=1)
    return (H1-L1)@n.T/C*1e3
for name,prob in (("row-major", raw.reshape(-1)),
                  ("col-major", raw.T.reshape(-1))):
    for oname,hp in (("nested",hp_n),("ring",hp_r)):
        dt=delay_for(hp)
        b=np.arange(-10.25,10.25,0.5)
        hist,_=np.histogram(dt,bins=b,weights=prob)
        hh=hist/hist.sum()
        peak=hh.max()
        # sharpness: fraction of mass in best 1-ms window (2 bins)
        best=max(hh[i]+hh[i+1] for i in range(len(hh)-1))
        print(f"{name:10s} {oname:7s} peak_bin={peak:.4f} best1ms={best:.4f} argmax_delay={b[np.argmax(hh)]+0.25:+.2f} ms")
