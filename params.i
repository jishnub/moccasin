
integer nr, smin, smax, lmin, lmax, ns, ordmax, ntot, nt
integer ordmin, nsplits, nord, nordmax, dsig
real*8 dt, freqmin, freqmax, sigmin, sigmax , offresonance 
character*80 mdidir, hmidir, freqmdidir, freqhmidir,QDPpath
character*80 outmdidir, outhmidir
logical FLOW_ANALYSIS, compute_norms
real*8 rotEarth, trackrate
parameter(rotEarth = 0.0317, trackrate = 0.453)!\mu Hz
real*8 pi, rsun, msun, Gcon
parameter (pi = acos(-1.0d0), Gcon = 6.678e-8)
character*3 instrument
parameter ( smin = 1, smax = 50, lmin = 1, lmax = 249)
parameter (nsplits = 36)
parameter (hmidir = '/scratch/jb6888/HMI')
parameter (mdidir = '/scratch/jb6888/MDI')
parameter (outhmidir = '/scratch/jb6888/HMIout')
parameter (outmdidir = '/scratch/jb6888/MDIout')
parameter (freqhmidir = '/scratch/jb6888/freqs')
parameter (freqmdidir = '/scratch/jb6888/freqs')
parameter (QDPpath = '/scratch/jb6888/QDP')

parameter(sigmin = 0.0, sigmax = 20.0, offresonance = 3.0, dsig = 1, freqmin = 1500.0, freqmax = 4000.0)

