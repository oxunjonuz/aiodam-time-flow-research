import numpy as np, math
RA=math.radians(134.79545455); DEC=math.radians(-69.79389373)
n=np.array([math.cos(DEC)*math.cos(RA), math.cos(DEC)*math.sin(RA), math.sin(DEC)])
C=299792458.0
VERT={"H1":np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6]),
      "L1":np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])}
d={k: float(v@n)/C for k,v in VERT.items()}
print("r.n/c  H1 %.6f ms  L1 %.6f ms"%(d["H1"]*1e3, d["L1"]*1e3))
print("H1 - L1 = %.4f ms"%( (d["H1"]-d["L1"])*1e3))
print("L1 - H1 = %.4f ms"%( (d["L1"]-d["H1"])*1e3))
# L1 arrived first => arrival time smaller at L1
