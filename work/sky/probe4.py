import numpy as np, math
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy import units as u

LAL = {  # official LALDetectors.h
 "H1": dict(lon=-2.08405676917, lat=0.81079526383, xaz=5.65487724844, yaz=4.08408092164,
            x=np.array([-0.22389266154,0.79983062746,0.55690487831]),
            y=np.array([-0.91397818574,0.02609403989,-0.40492342125])),
 "L1": dict(lon=-1.58430937078, lat=0.53342313506, xaz=4.40317772346, yaz=2.83238139666,
            x=np.array([-0.95457412153,-0.14158077340,-0.26218911324]),
            y=np.array([0.29774156894,-0.48791033647,-0.82054461286])),
}
def arm_vecs(lon, lat, yangle, xangle, yalt=0.0, xalt=0.0):
    resp = np.array([[-1,0,0],[0,0,0],[0,0,0]])
    rm2 = rotation_matrix(-lon*u.rad, 'z')
    rm1 = rotation_matrix(-1.0*(np.pi/2.0-lat)*u.rad, 'y')
    vecs = []
    for angle, azi in [(yangle, yalt), (xangle, xalt)]:
        rm0 = rotation_matrix(angle*u.rad, 'z')
        rmN = rotation_matrix(-azi*u.rad, 'y')
        rm = rm2 @ rm1 @ rm0 @ rmN
        vecs.append(rm @ np.array([-1,0,0]))
    return vecs  # [y, x]
for det,d in LAL.items():
    yv, xv = arm_vecs(d["lon"], d["lat"], d["yaz"], d["xaz"])
    print(det)
    print("  computed x:", np.round(xv,9), " LAL x:", d["x"], " diff", np.abs(xv-d["x"]).max())
    print("  computed y:", np.round(yv,9), " LAL y:", d["y"], " diff", np.abs(yv-d["y"]).max())
