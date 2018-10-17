import sys,os,subprocess
iterno=int(sys.argv[1]);nprocs=80;
llow=int(70);lhigh=int(149);nyears=int(4);styear=int(1);dell=int(0);instrum='HMI'
yearlen=int(4)
proc_no=int(os.environ['PBS_VNODENUM'])
jobno=(iterno-1)*nprocs+proc_no
def get_params(jobno):
  jobno_counter=0
  for year in xrange(styear,nyears+1,yearlen):
    for l in xrange(llow,lhigh+1):
        lprime = l + dell
        if jobno_counter==jobno: return year,l,lprime
        jobno_counter+=1
year,l,lprime = get_params(jobno)
os.chdir('/home/shravan/poloidal')
if not os.path.exists(os.path.join('/scratch/shravan',instrum,'processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))):
    print 'quitting, path not found:', os.path.join('/scratch/shravan',instrum,'processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))
    quit()
subprocess.call(["/home/shravan/poloidal/analyze",'--yearnum',str(year),'--ell',str(l),'--ellp', str(lprime),'--instrument',instrum,'--nyears',str(yearlen)])
