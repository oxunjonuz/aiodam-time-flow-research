import numpy as np, math
# Official LAL arm direction unit vectors, Earth-centered frame (LALDetectors.h)
ARM = {
 "H1": dict(x=np.array([-0.22389266154, 0.79983062746, 0.55690487831]),
            y=np.array([-0.91397818574, 0.02609403989, -0.40492342125])),
 "L1": dict(x=np.array([-0.95457412153, -0.14158077340, -0.26218911324]),
            y=np.array([ 0.29774156894, -0.48791033647, -0.82054461286])),
}
for k,v in ARM.items():
    print(k, "x norm", np.linalg.norm(v["x"]), "y norm", np.linalg.norm(v["y"]),
          "dot", float(v["x"]@v["y"]))

def gmst(gps):
    # IAU 1982 GMST (same family as LAL's gmst_accurate)
    # days since J2000.0 (2000-01-01 12:00:00 UTC = GPS 630763213 ... use unix)
    # GPS epoch 1980-01-06 00:00:00 UTC = unix 315964800
    unix = 315964800 + gps - 18.0   # GPS-UTC leap seconds = 18 at 2015
    jd = unix/86400.0 + 2440587.5
    T = (jd - 2451545.0)/36525.0
    g = 280.46061837 + 360.98564736629*(jd-2451545.0) + 0.000387933*T*T - T*T*T/38710000.0
    return math.radians(g % 360.0)

def fpfx(det, ra, dec, psi, t_gps):
    gha = gmst(t_gps) - ra
    cgha, sgha = math.cos(gha), math.sin(gha)
    cdec, sdec = math.cos(dec), math.sin(dec)
    cpsi, spsi = math.cos(psi), math.sin(psi)
    u = ARM[det]["x"]; v = ARM[det]["y"]
    D = 0.5*(np.outer(v,v) - np.outer(u,u))
    x = np.array([-cpsi*sgha - spsi*cgha*sdec,
                  -cpsi*cgha + spsi*sgha*sdec,
                   spsi*cdec])
    y = np.array([ spsi*sgha - cpsi*cgha*sdec,
                   spsi*cgha + cpsi*sgha*sdec,
                   cpsi*cdec])
    Fp = float(x@D@x - y@D@y)
    Fc = float(x@D@y + y@D@x)
    return Fp, Fc

RA = math.radians(134.79545455); DEC = math.radians(-69.79389373)
GPS = 1126259462.4
print()
for det in ("H1","L1"):
    for psi_deg in (0, 45, 90):
        Fp, Fc = fpfx(det, RA, DEC, math.radians(psi_deg), GPS)
        print(f"{det} psi={psi_deg:3d}  F+={Fp:+.4f} Fx={Fc:+.4f}  F+^2+Fx^2={Fp*Fp+Fc*Fc:.4f}")
