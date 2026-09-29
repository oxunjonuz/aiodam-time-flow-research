import time, math, numpy as np
from phenomxpy import IMRPhenomT
MC=28.096; Q=29/36; ETA=Q/(1+Q)**2; MTOT=MC/ETA**0.6
def build(mc, inc_deg, spins=(0.0,0.0)):
    eta = (mc/(mc/ETA**0.6))**1.0  # keep eta fixed, vary total mass from mc
    mtot = mc/ETA**0.6
    w=IMRPhenomT(eta=ETA,s1=[0,0,spins[0]],s2=[0,0,spins[1]],f_min=20.0,f_ref=20.0,
                 total_mass=mtot,distance=410.0,inclination=math.radians(inc_deg),
                 phi_ref=0.0,f_max=1024.0,delta_t=0.5/1024.0,delta_f=0.25,condition=True)
    fd=w.compute_fd_polarizations()
    return np.asarray(fd[0],dtype=complex), np.asarray(fd[1],dtype=complex)
t0=time.time(); hp,hc=build(MC,0.0); t1=time.time()
print("first call %.3f s, len %d"%(t1-t0,len(hp)))
t0=time.time()
for i in range(5): build(MC*(1+0.01*i),0.0)
print("5 more calls: %.3f s total, %.3f s each"%((time.time()-t0),(time.time()-t0)/5))
print("hp[100]",hp[100])
for inc in (0.0,30.0,60.0,90.0):
    a,b=build(MC,inc)
    print(f"inc={inc:5.1f} |hp|/|hp0|={np.abs(a[100])/np.abs(hp[100]):.6f} (1+cos^2)/2={(1+math.cos(math.radians(inc))**2)/2:.6f} |hc|/|hp0|={np.abs(b[100])/np.abs(hp[100]):.6f}")
