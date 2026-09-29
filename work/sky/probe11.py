import numpy as np, math, healpy as hp
C=299792458.0
H1=np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6])
L1=np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])
def gmst(gps):
    unix = 315964800 + gps - 18.0
    jd = unix/86400.0 + 2440587.5
    T=(jd-2451545.0)/36525.0
    g=280.46061837+360.98564736629*(jd-2451545.0)+0.000387933*T*T-T*T*T/38710000.0
    return math.radians(g%360.0)
GPS=1126259462.4
g=gmst(GPS)
print("GMST deg %.4f"%(math.degrees(g)))
m=hp.read_map('LALInference_skymap.fits',field=0,verbose=False)
npix=len(m); pix=np.arange(npix)
th,ph=hp.pix2ang(512,pix,nest=True)
ra=ph; dec=np.pi/2-th
ncel=np.stack([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)],axis=1)
# rotate celestial -> earth-fixed by GMST about z
cg,sg=math.cos(g),math.sin(g)
R=np.array([[cg,sg,0],[-sg,cg,0],[0,0,1]])
near=(R@ncel.T).T
proj=(H1-L1)@near.T/C*1e3
k=int(np.argmax(m))
print("ML pixel delay (H1-L1) %.4f ms"%proj[k])
b=np.arange(-10.05,10.05,0.2)
hist,_=np.histogram(proj,bins=b,weights=m)
hh=hist/hist.sum(); i=np.argmax(hh)
print("peak bin %.2f ms mass %.4f"%(b[i]+0.1,hh[i]))
for j in np.argsort(hh)[::-1][:5]: print("   %+.2f : %.4f"%(b[j]+0.1,hh[j]))
