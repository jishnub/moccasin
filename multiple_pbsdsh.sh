#!/bin/bash
#PBS -l nodes=5:ppn=2
#PBS -o  output-multipbsdsh1
#PBS -e  error-multipbsdsh1
#PBS -l walltime=48:00:00
cd $PBS_O_WORKDIR

instrument='HMI';
llow=1;lhigh=10;nyears=6;yearlength=6;styear=1;dell=1;
##lprimehigh=300;
#nparams=`python -c "print (-1 + $llow - $lhigh)*(-2 + $llow + $lhigh - 2*$lhigh)/2*$nyears"`
nparams=10 #`python -c "print ($lhigh-$llow+1)*($nyears//$yearlength)"`  
nprocs=10 ## nprocs * ppn
nloops=`python -c "import math;print int(math.ceil($nparams*1./$nprocs))"`
jobscript="$PBS_O_WORKDIR/jobscript$styear.py"
##jobscript="/home/shravan/poloidal/jobscript$styear.py"
execpath="/home/shravan/poloidal/analyze"

echo "import sys,os,subprocess" > $jobscript
echo "iterno=int(sys.argv[1]);nprocs=$nprocs;" >> $jobscript
echo "llow=int($llow);lhigh=int($lhigh);nyears=int($nyears);styear=int($styear);dell=int($dell);instrum='$instrument'" >> $jobscript ##lprimehigh=int($lprimehigh);
echo "yearlen=int($yearlength)" >> $jobscript ##lprimehigh=int($lprimehigh);
echo "proc_no=int(os.environ['PBS_VNODENUM'])" >> $jobscript
echo "jobno=(iterno-1)*nprocs+proc_no" >> $jobscript
##echo "print jobno, proc_no" >>$jobscript
#echo "print iterno,nprocs,proc_no,jobno" >> $jobscript
echo "def get_params(jobno):" >> $jobscript
echo "  jobno_counter=0" >> $jobscript
echo "  for year in xrange(styear,nyears+1,yearlen):" >> $jobscript
echo "    for l in xrange(llow,lhigh+1):" >> $jobscript
##echo "      for lprime in xrange(l+dell+1,l+dell+1):" >> $jobscript #lhigh+1
##echo "      lprime =l+dell" >> $jobscript #lhigh+1
##echo "      for lprime in xrange(l,lhigh+1):" >> $jobscript #lhigh+1
##echo "        print jobno,jobno_counter,l,lprime,year" >> $jobscript
echo "        lprime = l + dell" >> $jobscript
echo "        if jobno_counter==jobno: return year,l,lprime" >> $jobscript
echo "        jobno_counter+=1" >> $jobscript
echo "year,l,lprime = get_params(jobno)" >> $jobscript
echo "os.chdir('$PBS_O_WORKDIR')" >> $jobscript
##echo os.path.join('{0:03d}_{1:03d}_{0:02d}'.format(l,lprime,year))
#echo "if not os.path.exists(os.path.join('/scratch/shravan',\"$instrument\",'processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))):" >> $jobscript
#echo "if not os.path.exists(os.path.join('/scratch/shravan',\"${instrument}\",'processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))):" >> $jobscript
#echo "print os.path.join('/scratch/shravan/HMI/processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))"
echo "if not os.path.exists(os.path.join('/scratch/shravan',instrum,'processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))):" >> $jobscript
##echo "if not os.path.exists(os.path.join('/scratch/shravan/MDI/processed','{0:02d}_{1:03d}'.format(year,l))):" >> $jobscript
##echo "if not os.path.exists(os.path.join('/scratch/shravan/HMI/processed','{0:02d}_{1:03d}'.format(year,l))):" >> $jobscript
echo "    print 'quitting, path not found:', os.path.join('/scratch/shravan',instrum,'processed','{0:03d}_{1:03d}_{2:02d}'.format(l,lprime,year))" >> $jobscript
echo "    quit()" >> $jobscript
##echo "subprocess.call([\"${execpath}\ --ell",str(l),"--yearnum",str(year),"--ellp",str(lprime),"nyears",str($yearlength),"--instrument",$instrument])" >> $jobscript
##echo 'print (["/home/shravan/poloidal/analyze","--yearnum",str(year),"--ell",str(l),"--ellp", str(lprime),"--instrument",instrum,"--nyears",str(yearlen)])' >> $jobscript
##echo 'subprocess.call(["/home/shravan/poloidal/analyze","--yearnum",str(year),"--ell",str(l),"--ellp", str(lprime),"--instrument",instrum,"--nyears",str(yearlen)])' >> $jobscript

echo "subprocess.call([\"${execpath}\",'--yearnum',str(year),'--ell',str(l),'--ellp', str(lprime),'--instrument',instrum,'--nyears',str(yearlen)])" >> $jobscript
##echo "subprocess.call([\"${execpath}\",'--yearnum',str(year),'--ell',str(l),'--ellp', str(lprime),'--instrument',instrum,'--nyears',str(yearlen),'--compute_norms'])" >> $jobscript

for iterno in `seq 1 $nloops`
do
	/usr/local/bin/pbsdsh /home/shravan/anaconda2/bin/python $jobscript $iterno
done
#rm $jobscript
