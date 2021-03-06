
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
parameter ( smin = 1, smax = 7, lmin = 1, lmax = 10)
parameter (nsplits = 36)
parameter (hmidir = '/scratch/jishnu/data/HMI/data')
parameter (mdidir = '/scratch/jishnu/data/MDI/data')
parameter (outhmidir = '/scratch/shravan/HMI')
parameter (outmdidir = '/scratch/shravan/MDI')
parameter (freqhmidir = '/scratch/shravan/freqs')
parameter (freqmdidir = '/scratch/shravan/freqs')

parameter(sigmin = 0, sigmax = 200.0, offresonance = 6.0, dsig = 1, freqmin = 1400.0, freqmax = 4600.0)

