from astropy.io import fits
import numpy as np; import array
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy
import math
from numpy.linalg import inv
import numpy.matlib

def dfac(n):
    a = np.zeros_like(n)
    for j in range(0,n.size):
      a[j] = np.prod(np.arange(n[j],1,-2))
    return a

def fac(n):
    a = np.zeros_like(n)
    for j in range(0,n.size):
      a[j] = math.factorial(n[j])
    return a


elmin = 1
elmax = 10
dellmin = 1
dellmax =1
instrument = "HMI"
track = 430
smin =0
smax = 7
perform_inversion = True
sigma_min = 1 # micro-Hertz
sigma_max = 200 # micro-Hertz
nyears = 6 # each "year" is 360 days
startyear = 1

r0 = 0.15
deltar = 0.05
reg = 1e6
dnu = 1e6/(nyears*360.0*86400.0) # in mu Hz
freqfloor = int(np.floor(sigma_min/dnu) +1)
r=np.squeeze(fits.open("radius.fits")[0].data)
trackch =str("{:03d}".format(track));
lminst =str("{:03d}".format(elmin));
lmaxst =str("{:03d}".format(elmax));
direct0 = '/scratch/shravan/HMI'
direct = direct0 + '/tracking' + trackch
ns = (smax+1)**2
dellst =str("{:01d}".format(dellmin));
kern = np.squeeze(fits.open(instrument+"_kernels_"+lminst +"_to_"+lmaxst+"_dell_"+dellst+".fits")[0].data)
indices = np.loadtxt(instrument+"_indices_"+lminst +"_to_"+lmaxst+"_dell_"+dellst)

for dell in range(dellmin+1, dellmax+1):
    dellst =str("{:01d}".format(dell));
    kerntemp = np.squeeze(fits.open(instrument+"_kernels_"+lminst +"_to_"+lmaxst+"_dell_"+dellst+".fits")[0].data)
    g = np.loadtxt(instrument+"_indices_"+lminst +"_to_"+lmaxst+"_dell_"+dellst)
    indices = np.r_[indices, g]
    kern = np.r_[kern, kerntemp]

h = kern.shape
nkerns = h[0]

dr = np.zeros_like(r)

nr = r.shape
nr = nr[0]

dr[0:nr-2] = r[1:nr-1] - r[0:nr-2]
dr[nr-1] = dr[nr-2]

target = np.exp(-(r-r0)**2.0/(2.0*deltar**2))
target = target/np.sum(target*dr)

A = np.zeros((nkerns,nkerns))
rhs = np.zeros((nkerns,1))


for j in range(0,nkerns):
    temp = dr*kern[j,:]
    rhs[j] = np.sum(temp*target)
    for i in range(j,nkerns):
        A[i,j] = np.sum(temp*kern[i,:])
        A[j,i] = A[i,j]


parity = np.zeros((nkerns,1))
coeffstor = np.zeros((nkerns,smax+1))
coeffspol = np.zeros((nkerns,smax+1))

elldiff =  indices[:,2] - indices[:,0]


for s in range(smin,smax+1):
     
     parity = np.mod(elldiff + s,2)
     hh = np.where(abs(elldiff) <= s)[0]
     if (hh.size == 0):
          continue     

     tind = hh[np.where(parity[hh] == 1)[0]]
     #tind = tind[hh]

     pind = hh[np.where(parity[hh] == 0)[0]]
     #pind = np.where(parity == 0 and abs(elldiff) <= s)

     sumdiff1 = s + elldiff
     sumdiff2 = s - elldiff

     if (tind.size >0):
         factor = (1-2*np.mod((sumdiff1[tind] - 1)/2,2)) * dfac(sumdiff1[tind]) * dfac(sumdiff2[tind])/(fac(sumdiff1[tind])*fac(sumdiff2[tind]))**0.5
         rhstor = rhs[tind,0] * factor
         Ator = A[np.ix_(tind, tind)]
         coeffstor[tind,s] =  np.matmul(inv(Ator + reg * np.eye(tind.size)), rhstor)
#         if (s==2):
 #           func = np.matmul(np.squeeze(coeffstor[tind,s]),kern[tind,:])
  #          plt.plot(r,target/target.max()); plt.plot(r,func/func.max()); plt.show()
   #         stop

     if (pind.size >0):
         facpol = (1-2*np.mod(sumdiff1[pind]/2,2)) * elldiff[pind] * dfac(sumdiff1[pind]-1) * dfac(sumdiff2[pind]-1)/(fac(sumdiff1[pind])*fac(sumdiff2[pind]))**0.5
         rhspol = rhs[pind,0] * facpol
         Apol = A[np.ix_(pind, pind)]
         coeffspol[pind,s] =  np.matmul(inv(Apol + reg * np.eye(pind.size)), rhspol)
         
if (perform_inversion):
   nfreq = int(np.floor((sigma_max - sigma_min)/dnu)) + 2
   a = np.zeros([nfreq,60,ns],'complex'); powpos = 0; powneg = 0;
   nus  = (np.arange(nfreq) + freqfloor)*dnu 
   noitoroidal = np.zeros([nfreq,ns])
   noipoloidal = np.zeros([nfreq,ns])
   toroidal = np.zeros([nfreq,ns], dtype = complex)
   poloidal = np.zeros([nfreq,ns], dtype = complex)
   nors = np.zeros([60]);nord = np.zeros([60]); nordp = np.zeros([60])
   stry1 = str("{:02d}".format(5*(startyear-1)+1));
   stry2 = str("{:02d}".format(5*(startyear+nyears-1)));
   stryear = '_year_'+stry1+'_'+stry2

   for dell in range(dellmin, dellmax+1):
      for ell in range(elmin, elmax+1-dell):
        print "ell:", ell
        ellc =str("{:03d}".format(ell))
        ellp = ell + dell
        ellpc =str("{:03d}".format(ellp))
        allind = np.where(indices[:,0] == ell)[0]
        allind = allind[np.where(indices[allind,2] == ellp)]
        
        te = fits.open(direct+'/bcoef_l_'+ellc +'_lp_'+ellpc+stryear+'.fits')[0].data
        noit = fits.open(direct+'/noise_l_'+ellc +'_lp_'+ellpc+stryear+'.fits')[0].data
        nfrequ = te.shape[1]
        te = te[0,:,:,:]+1j*te[1,:,:,:]
        f = open(direct+'/frequency_metadata_l_'+ellc +'_lp_'+ellpc+stryear, 'r')
        j=-1
        k= -1

        for line in f:

            j = j+1
            if (j==0):
                line = line.strip()
                columns = line.split()
                dnu = float(columns[0])
            if (j <= 3):
                continue 
            k = k+1
            line = line.strip()
            columns = line.split()
            freqdiff = np.float(columns[5]) - np.float(columns[2]) 
            if (np.absolute(freqdiff) < sigma_min or  np.absolute(freqdiff) > sigma_max):
                 continue

            fst = int(np.floor(np.absolute(freqdiff)/dnu)) - freqfloor 
            fend = np.minimum(fst + nfrequ, nfreq)
            nord = int(columns[1])
            nordp = int(columns[4])
            nind = np.where(indices[allind,1] == nord)[0]
            freql = nfrequ + (fend-fst - nfrequ)
            #print fst,fend,nus[fst],nus[fend-1],fend-fst,nfrequ,freql,nus.shape

#            if (fend > nfreq-1):
 #               print "Frequency range too high, skipping, ell =", ell, "dell =", dell, "n = ", nord, "n' = ", nordp, "freq. difference = ", freqdiff
  #              continue

            coefind = allind[np.where(indices[allind[nind],3] == nordp)][0]
            for s in range(smin,smax+1):
                for t in range(-s,s+1):
                    ind = s**2 + s + t
                   
                    poloidal[fst:fend,ind] = poloidal[fst:fend,ind] +  te[0:freql,k,ind]*coeffspol[coefind,s]*(1.0 + 0.0j)
                    noipoloidal[fst:fend,ind] = noipoloidal[fst:fend,ind] +  noit[0:freql,k,ind]*coeffspol[coefind,s]

                    toroidal[fst:fend,ind] = toroidal[fst:fend,ind] +  te[0:freql,k,ind]*coeffstor[coefind,s]*(1.0 + 0.0j)
                    noitoroidal[fst:fend,ind] = noitoroidal[fst:fend,ind] +  noit[0:freql,k,ind]*coeffstor[coefind,s]

        f.close()

