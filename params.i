
integer nr, smin, smax, lmin, lmax, ns, ordmax, ntot, nt
integer ordmin, nsplits, nord, sigmin, sigmax, nordmax
real*8 dt, freqmin, freqmax !, domega
logical MODE_ANALYSIS
real*8 rotEarth, trackrate
parameter(rotEarth = 0.0317, trackrate = 0.440)!\mu Hz
real*8 pi, rsun, msun, Gcon
parameter (pi = acos(-1.0d0), Gcon = 6.678e-8)
character*3 instrument
parameter(instrument = 'HMI')
parameter ( smin = 4, smax = 30, lmin = 44, lmax = 255)
parameter ( nsplits = 36)
!, nord = ordmax-ordmin+1, nordmax = 7) 
parameter(sigmin = 1, sigmax = 75, freqmin = 1400.0, freqmax = 4600.0)
!nt = 1036800,int(8.64*7.2*0.25),int(8.64*7.2*6)
!character*(500) modeparams
!parameter(modeparams='mdi.vw_V_sht_modes.19960501_000000_TAI.0.300.518400.m10qr.1216.36')

