import numpy as np, math, healpy as hp
from astropy.io import fits
C=299792458.0
H1=np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6])
L1=np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])
def gmst(gps):
    unix=315964800+gps-18.0; jd=unix/86400.0+2440587.5; T=(jd-2451545.0)/36525.0
    g=280.46061837+360.98564736629*(jd-2451545.0)+0.000387933*T*T-T*T*T/38710000.0
    return math.radians(g%360.0)
G=gmst(1126259462.4)
raw=np.asarray(fits.open('LALInference_skymap.fits')[1].data['PROB']).astype(float)
npix=3145728; pix=np.arange(npix)
th,ph=hp.pix2ang(512,pix,nest=True)
ra=ph; dec=np.pi/2-th
ncel=np.stack([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)],axis=1)
def rot(n,g):
    c,s=math.cos(g),math.sin(g); R=np.array([[c,s,0],[-s,c,0],[0,0,1]]); return (R@n.T).T
cands={'row':raw.reshape(-1),'col':raw.T.reshape(-1)}
b=np.arange(-10.05,10.05,0.1)
for cn,prob in cands.items():
    for rn,g in (('none',0.0),('+G',G),('-G',-G)):
        near = ncel if g==0.0 else rot(ncel,g)
        proj=(H1-L1)@near.T/C*1e3
        h,_=np.histogram(proj,bins=b,weights=prob); hh=h/h.sum()
        i=np.argmax(hh)
        k=int(np.argmax(prob))
        print(f'{cn:4s} {rn:5s} peak {hh[i]:.4f}@{b[i]+0.05:+.2f}ms   MLdelay {proj[k]:+.2f}ms  ML ra {np.degrees(ra[k]):.1f} dec {np.degrees(dec[k]):.1f}')
