import numpy as np, math, time
# LAL geometry (LALDetectors.h, git.ligo.org/lscsoft/lalsuite)
DET = {
 "H1": dict(lon=-2.08405676917, lat=0.81079526383, elev=142.554,
            xaz=5.65487724844, yaz=4.08408092164),
 "L1": dict(lon=-1.58430937078, lat=0.53342313506, elev=-6.574,
            xaz=4.40317772346, yaz=2.83238139666),
}
def dircos(lon, lat, az):
    # unit vector pointing along an arm, Earth-centered frame
    # from LAL: x = cos(lat)cos(lon-az)... use standard
    return np.array([math.cos(lat)*math.cos(lon)*math.cos(az)+math.sin(lon)*math.sin(az),
                     math.cos(lat)*math.sin(lon)*math.cos(az)-math.cos(lon)*math.sin(az),
                     math.sin(lat)*math.cos(az)])
# Use LAL's published arm direction vectors to validate
LAL_X = {"H1": np.array([-0.22389266154, 0.79983062746, 0.55690487831]),
         "L1": None}
print("H1 x-arm from formula:", dircos(DET["H1"]["lon"],DET["H1"]["lat"],DET["H1"]["xaz"]))
print("H1 x-arm from LAL     :", LAL_X["H1"])
