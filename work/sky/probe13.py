import numpy as np, math, h5py, os
from scipy.signal import butter, sosfiltfilt
from phenomxpy import IMRPhenomT
FS=16384; GPS0=1126259447; T_EVENT=1126259462.4-GPS0
FMIN,FMAX=20.0,300.0; MC=28.096; Q=29/36; ETA=Q/(1+Q)**2; MTOT=MC/ETA**0.6
DF=0.25; FMAX_FD=1024.0
DATA="/work/time_flow_research_20260927/data"
FILES={"H1":"H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5","L1":"L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}
ARM={"H1":(np.array([-0.22389266154,0.79983062746,0.55690487831]),
           np.array([-0.91397818574,0.02609403989,-0.40492342125])),
     "L1":(np.array([-0.95457412153,-0.14158077340,-0.26218911324]),
           np.array([ 0.29774156894,-0.48791033647,-0.82054461286]))}
VERT={"H1":np.array([-2.16141492636e6,-3.83469517889e6,4.60035022664e6]),
      "L1":np.array([-7.42760447238e4,-5.49628371971e6,3.22425701744e6])}
C=299792458.0
def gmst(gps):
    unix=315964800+gps-18.0; jd=unix/86400.0+2440587.5; T=(jd-2451545.0)/36525.0
    return math.radians((280.46061837+360.98564736629*(jd-2451545.0)+0.000387933*T*T-T*T*T/38710000.0)%360.0)
def fpfx(det,ra,dec,psi,gps):
    gha=gmst(gps)-ra
    cgha,sgha=math.cos(gha),math.sin(gha); cdec,sdec=math.cos(dec),math.sin(dec)
    cpsi,spsi=math.cos(psi),math.sin(psi)
    x,y=ARM[det]; D=0.5*(np.outer(y,y)-np.outer(x,x))
    v=np.array([-cpsi*sgha-spsi*cgha*sdec, -cpsi*cgha+spsi*sgha*sdec, spsi*cdec])
    w=np.array([ spsi*sgha-cpsi*cgha*sdec,  spsi*cgha+cpsi*sgha*sdec, cpsi*cdec])
    return float(v@D@v-w@D@w), float(v@D@w+w@D@v)
def load(det):
    with h5py.File(os.path.join(DATA,FILES[det]),"r") as h:
        return h["strain/Strain"][:], int(h["meta/GPSstart"][()])
def hp_filter(x,fs,f0=15.0): return sosfiltfilt(butter(4,f0/(fs/2.0),btype="highpass",output="sos"),x)
def welch(x,fs,nper):
    step=int(nper*0.5); w=np.hanning(nper); norm=2.0/(fs*np.sum(w**2))
    segs=np.array([np.abs(np.fft.rfft(x[i:i+nper]*w))**2*norm for i in range(0,len(x)-nper,step)])
    return np.fft.rfftfreq(nper,1.0/fs), np.median(segs,axis=0)/math.log(2.0)
def psd_on(f,fp,Pp): return np.exp(np.interp(np.log(np.clip(f,1e-3,None)),np.log(np.clip(fp,1e-3,None)),np.log(Pp)))
def model(mc,inc):
    mtot=mc/ETA**0.6
    w=IMRPhenomT(eta=ETA,s1=[0,0,0],s2=[0,0,0],f_min=FMIN,f_ref=FMIN,total_mass=mtot,
                 distance=410.0,inclination=math.radians(inc),phi_ref=0.0,f_max=FMAX_FD,
                 delta_t=0.5/FMAX_FD,delta_f=DF,condition=True)
    fd=w.compute_fd_polarizations()
    return np.asarray(fd[0],dtype=complex), np.asarray(fd[1],dtype=complex)
seg_n=4*FS; i0=int(round((T_EVENT-2.0)*FS))
freqs=np.fft.rfftfreq(seg_n,1.0/FS)
m=(freqs>=FMIN)&(freqs<=FMAX)
hp,hc=model(MC,0.0)
f_model=np.arange(len(hp))*DF
Hp=np.interp(freqs,f_model,hp.real)+1j*np.interp(freqs,f_model,hp.imag)
Hc=np.interp(freqs,f_model,hc.real)+1j*np.interp(freqs,f_model,hc.imag)
RA=math.radians(134.79545455); DEC=math.radians(-69.79389373)
print("antenna responses at ML sky position:")
for psi_deg in (0,30,60,90,120,150):
    psi=math.radians(psi_deg)
    r={}
    for det in ("H1","L1"):
        Fp,Fc=fpfx(det,RA,DEC,psi,1126259462.4)
        # projected template
        Ht=Fp*Hp+Fc*Hc
        Sn=psd_on(freqs,*welch(hp_filter(np.concatenate([load(det)[0][:8*FS],load(det)[0][-8*FS:]]),FS),FS,FS))
        rho=math.sqrt(4.0*DF*float(np.sum(np.abs(Ht[m])**2/Sn[m])))
        r[det]=rho
    print(f"psi={psi_deg:3d}  rho_opt H1={r['H1']:7.3f} L1={r['L1']:7.3f}  ratio={r['H1']/r['L1']:.4f}")
print("published ratio 19.5/13.3 =", 19.5/13.3)
