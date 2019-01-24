
integer nr, smin, smax, lmin, lmax, ns, ordmax, ntot, nt
integer ordmin, nsplits, nord, nordmax, dsig
real*8 dt, freqmin, freqmax, sigmin, sigmax , offresonance 
character*80 mdidir, hmidir, freqmdidir, freqhmidir
character*80 outmdidir, outhmidir
logical FLOW_ANALYSIS, compute_norms
real*8 rotEarth, trackrate
parameter(rotEarth = 0.0317, trackrate = 0.430)!\mu Hz
real*8 pi, rsun, msun, Gcon
parameter (pi = acos(-1.0d0), Gcon = 6.678e-8)
character*3 instrument
parameter ( smin = 1, smax = 1, lmin = 1, lmax = 255)
parameter (nsplits = 36)!, nr = 7305)
parameter (hmidir = '/scratch/jb6888/HMI')
parameter (mdidir = '/scratch/jb6888/MDI')
parameter (outhmidir = '/scratch/jb6888/HMIout')
parameter (outmdidir = '/scratch/jb6888/MDI')
parameter (freqhmidir = '/scratch/jb6888/freqs')
parameter (freqmdidir = '/scratch/jb6888/freqs')

parameter(sigmin = 0, sigmax = 200.0, offresonance = 6.0, dsig = 1, freqmin = 1400.0, freqmax = 4600.0)

