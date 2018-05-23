MODULE DATA_ANALYSIS

! use kernels
 implicit none
 include 'params.i'
 include 'fftw3.f'


 real*8 dnu!, rotEarth
 !parameter(rotEarth = 0.0317)!\mu Hz
 real*8, allocatable, dimension(:,:) :: nu
 real*8, allocatable, dimension(:,:,:) :: rleaks, horleaks
 real*8 freqnorm, rhonorm, speednorm !pi, rsun, msun, Gcon, 
 !parameter (pi = acos(-1.0d0), Gcon = 6.678e-8)

 complex*16, allocatable, dimension(:,:) :: shtdata, shtdatap

 type :: Cell
  complex*16, allocatable, dimension(:) :: tee
 end type  
 
 type(Cell), allocatable, dimension(:,:,:) :: b_coefs
! type(Cell), allocatable, dimension(:,:) :: corr, den

 type :: wigcell
  real*8, allocatable, dimension(:,:) :: sub3j
 end type  
 
 type(wigcell) :: wig(smin:smax)

 type :: Cell2
  real*8, allocatable, dimension(:) :: tee, em
 end type  

 type(Cell2), allocatable, dimension(:,:,:) :: noise, noisevar

 !type(Cell)
 !type :: Cell
 ! complex*16, allocatable, dimension(:,:) :: sub
 !end type  
 
 !type(Cell) :: Zmatrix(smin:smax, ordmin:ordmax), symbols(smin:smax)

 type :: subcell
  real*8, allocatable, dimension(:,:) :: mn
 end type 
 
 type(subcell) :: nus(lmin:lmax)

 type :: sub1
  real*8, allocatable, dimension(:) :: ords, ems
 end type 
  
 type(sub1) :: freqnu(lmin:lmax), en(lmin:lmax), a1(lmin:lmax),&
  amps(lmin:lmax), fwhm(lmin:lmax), backg(lmin:lmax), leakage(0:304)


 logical existence(lmin:lmax,0:24)
! type :: kern 
!  complex*16, allocatable, dimension(:,:) :: submat
! end type  
 
! type(kern) :: 



Contains
!------------------------------------------------------------------------




!------------------------------------------------------------------------
function Pf(l, N)

  implicit none
  integer i, j, twoellplus1, l, N
  real*8 P(-l:l,0:N), Lsq, em(-l:l), Pp(-l:l), c0
  real*8 mbyL(-l:l), PdP(-l:l), cj, Pf(-l:l,1:N) !Podd(-l:l,1:floor(N/2.0)+1), 

! Note N <= 2*l


  Lsq = (l*(l+1))**0.5
  twoellplus1 = (2*l+1)
  
  do i=-l,l
   em(i) = i
  enddo 

  P(:,0) = l
  P(:,1) = em
  mbyL = em/Lsq

  P(:,2) = Lsq**2.d0/(l-0.5d0) * legendre(2,mbyL,twoellplus1)

  do i=3,N
    Pdp = Lsq*legendre(i,mbyL,twoellplus1)
    
    c0 = sum(Pdp)/(l*(2*l+1))
    Pp = Pdp - c0 * P(:,0)

    do j=1,i-1
        cj = sum(Pdp*P(:,j))/sum(P(:,j)**2.)
        Pp = Pp - cj * P(:,j)
    enddo
    P(:,i) =  l * Pp/Pp(l)
  enddo

!Podd = P(:,1:N:2)

  Pf = P(:,1:N)

end function Pf

!------------------------------------------------------------------------

function legendre(ell, arg, larg)

  implicit none
  integer i, larg, ell
  real*8 arg(larg), Pmin2(larg), Pmin1(larg), legendre(larg)

  Pmin2 = 1
  Pmin1 = arg
  do i=2,ell
   legendre = 1.d0/dble(i) * ((2*i-1)*arg*Pmin1 - (i-1)*Pmin2)
   Pmin2 = Pmin1 
   Pmin1 = legendre
  enddo

end function legendre

!------------------------------------------------------------------------

subroutine construct_frequencies(inpfile,ll,llp,ynum)

  implicit none
  integer i, l, npoints, ind_data(1:25), j, ndata,llp
  integer minord, maxord, nrecs, ell, ellp, ll, ynum, outp(1)

  character*(*) inpfile
  character*3 ellc, ellpc
  character*2 yc
  !parameter(npoints = 2159)
  !real*8 splits(24,npoints)
  real*8, dimension(:), allocatable :: ells, ensall, freqs, amplitudes, &
                                       fwhmall, backgall!, a1, a3, a5
  !real*8 acoefs(0:2*nsplits,npoints)
  real*8, allocatable, dimension(:,:) :: Po, splits, acoefs, adata
  integer, allocatable, dimension(:) :: indicesn
  logical lexist

  l = ll
  write(ellc,'(I3.3)') ll
  write(ellpc,'(I3.3)') llp
  write(yc,'(I2.2)') ynum
  inquire(file='/home/shravan/QDP/lengs_'//ellc//'_'//ellpc//'_'//yc, exist=lexist)
  if (lexist) call system("rm /home/shravan/QDP/lengs_"//ellc//"_"//ellpc//'_'//yc)

  call system("awk '{ print NF }' "//inpfile//" > /home/shravan/QDP/lengs_"//ellc//"_"//ellpc//'_'//yc)

  open(52,file='/home/shravan/QDP/lengs_'//ellc//'_'//ellpc//'_'//yc,&
           status='old',position='rewind',action='read')
  read(52,*) nrecs
  close(52)

  call system("awk 'END {print NR}' "//inpfile//" > /home/shravan/QDP/lengs_"//ellc//"_"//ellpc//'_'//yc)

  open(52,file='/home/shravan/QDP/lengs_'//ellc//'_'//ellpc//'_'//yc,&
           status='old',position='rewind',action='read')
  read(52,*) npoints
  close(52)

  

  allocate(ells(npoints), ensall(npoints), freqs(npoints),amplitudes(npoints),&
        fwhmall(npoints), backgall(npoints), splits(nrecs,npoints), &
        acoefs(1:nsplits,npoints))

  open(334,file=inpfile,status='old',position='rewind',action='read')
  do i=1,npoints
   read(334, *) splits(:,i)
  enddo
  close(34)

  ells = splits(1,:)
  ensall = splits(2,:)
  freqs = splits(3,:)
  amplitudes = splits(4,:)
  fwhmall = splits(5,:)
  backgall = splits(6,:)
! x (7), 8, 9, 10, 11, sigma x (12)

  do i=1,nsplits
   acoefs(i,:) = splits(12+i,:) * 1e-3
  enddo
!  a1 = splits(13,:) * 1e-3 !a coeffs are read in at nHz. converting to \mu Hz
!  a3 = splits(15,:) * 1e-3
!  a5 = splits(17,:) * 1e-3

!%a6 (18), sigmas (6) - total 24

!  N = 3; !% maximum number of odd splits

  do l = ll-6, ll+ 6!lmin,lmax

    ndata = 0
    do i=1,npoints
     if ((int(ells(i)) == l)) then ! .and. (int(ensall(i)) .ne. 0) .and. (int(freqs(i)) .ne. 0)
         ndata = ndata + 1
         ind_data(ndata) = i
     endif
    enddo

    if (allocated(adata)) deallocate(adata)
    if (allocated(Po)) deallocate(Po)
    if (allocated(indicesn)) deallocate(indicesn)
    if (allocated(freqnu(l)%ords)) deallocate(freqnu(l)%ords)
    if (allocated(en(l)%ords)) deallocate(en(l)%ords)
    if (allocated(amps(l)%ords)) deallocate(amps(l)%ords)
    if (allocated(fwhm(l)%ords)) deallocate(fwhm(l)%ords)
    if (allocated(backg(l)%ords)) deallocate(backg(l)%ords)

    minord = int(ensall(ind_data(1)))
    maxord = int(ensall(ind_data(ndata)))

!     print *,minord,maxord,l, nsplits,ensall(ind_data(1)),ensall(ind_data(ndata)),ndata
    allocate(freqnu(l)%ords(minord:maxord), en(l)%ords(minord:maxord), amps(l)%ords(minord:maxord), &
      fwhm(l)%ords(minord:maxord), backg(l)%ords(minord:maxord), a1(l)%ords(minord:maxord),&
      adata(1:nsplits,minord:maxord),Po(-l:l, 1:Nsplits),indicesn(1:ndata))
!Po(-l:l, 1:(Nsplits/2+1))

    en(l)%ords(:) = -1.d0
    indicesn = int(ensall(ind_data(1:ndata)))
    en(l)%ords(indicesn) = ensall(ind_data(1:ndata))
    freqnu(l)%ords(indicesn) = freqs(ind_data(1:ndata))
    amps(l)%ords(indicesn) = amplitudes(ind_data(1:ndata))
    fwhm(l)%ords(indicesn) = fwhmall(ind_data(1:ndata))
    backg(l)%ords(indicesn) = backgall(ind_data(1:ndata))

    adata(:,indicesn) = acoefs(:,ind_data(1:ndata))
    a1(l)%ords(indicesn) = acoefs(1,ind_data(1:ndata))
    Po(-l:l,:) = Pf(l, Nsplits)

    if (l == ll)  then
     ordmin = int(en(l)%ords(minord - 1 +minloc(abs(freqmin - freqnu(l)%ords(:)),1)))
     ordmax = int(en(l)%ords(minord-1+minloc(abs(freqmax - freqnu(l)%ords(:)),1)))
    endif

    allocate(nus(l)%mn(-l:l,minord:maxord))

!    print *,minord,maxord, ordmin,ordmax
    do i=minord,maxord
     nus(l)%mn(:,i) = freqnu(l)%ords(i) 
     !do j=1,nsplits/2

     do j=-l,l ! SYNODIC CORRECTION TO THE ROTATION RATE - ONLY THE a1 COEFFICIENT IS ALTERED
       nus(l)%mn(j,i) = nus(l)%mn(j,i) + (adata(1,i) - rotEarth)* Po(j,1)!  Po(j,1)  - rotEarth
     enddo

     do j=2,nsplits ! NOW ADDING THE REST OF THE SPLITTING COEFFICIENTS
      nus(l)%mn(:,i) = nus(l)%mn(:,i) + adata(j,i) * Po(:,j)!a1_data(i)*Po(:,1) + a3_data(i)*Po(:,2) &
             !+ a5_data(i)*Po(:,3)
     enddo

    enddo
    deallocate(adata)

  enddo
  
!  open(132,file='frequencies.comparison',action='write',position='rewind',status='replace')
!  do i=-120,120
!   write(132,*) nus(120)%mn(i,:)
!  enddo
!  close(132)
!  print *,en(120)%ords
! stop
 
  deallocate(acoefs, ensall, freqs, amplitudes, fwhmall, backgall, Po)
end subroutine construct_frequencies

!------------------------------------------------------------------------

subroutine analyzethis (yearnum, ell, ellp, nyears)

 Use Ziggurat
 implicit none

 integer order, ell, slow, s, m, mp, maxoff, offseti, yearnum, eltp, i
 integer indm, indmp, indm_nu, indmp_nu, ind_nu, t, minoff,  sigind, seed
 integer ordlow, ordhigh,sigoff, siglow, sighigh, ellp, ellt, nyears
 integer emp, emdp, mmin, mmax, signemp, signempt, shif, dl, dm, dmp
 integer minoffp, maxoffp, eldp, delm, dell,lolim,hilim,refind,indst
 integer intarr(1:10000), j, k, ierr, emdpmin, emdpmax, reals, nreals
 parameter(delm = 15, dell = 6, nreals = 10)
 integer, allocatable, dimension(:) :: minoffarr, maxoffarr, seeds

 character*4 daynum
 character*3 lch, lch_p, trackch
 character*2 ynum, ynum2
! character*1 omin, omax

 real*8 hilim_mp, lolim_mp, lolim_m, hilim_m, gamnlp, cons, con0, defcon
 !real*8 hilim_empt, lolim_empt, lolim_emp, hilim_emp
 real *8 omega(nt) ,dnu, contwopi, sig, nulmn, nulmnt, gamnl, consp
 real*8 nus_m, nus_mp, tstart, tfin, test(-1:1), nj(-6:6), con,noisevar2(nreals)
 real*8  den(smin:smax), conr, mix(-6:6), fmode, compar, ref, temp1
 real*8 noisehalf(1:10000),noiseother(1:10000),noisevar1(nreals,smin:smax)
 complex*16, dimension(1:10000,1:nreals) :: temphalf, tempother
 real*8, dimension(1:40000) :: vec
 real*8, dimension(:), allocatable :: tempm
 real*8, allocatable, dimension(:,:) :: leaks, mask, tempr, leaksp, norms
 real*8, allocatable, dimension(:,:,:) :: abspow
 complex*16 conp
 complex*16, dimension(:,:,:), allocatable :: limitpow
 complex*16, allocatable, dimension(:,:) ::  prefac
 logical compute_noise, lexist, compute_ab_initio

 
! defcon = 1.6e-11 / nyears  
! if (instrument == 'HMI') defcon = 2.95e-11 / nyears  
 contwopi  = 1.d0/(2*pi)**2
 !compute_noise = .false.
 compute_noise = .true.
 compute_ab_initio = .true.
 if (yearnum == 1) compute_noise = .true.

 write(trackch,'(I3.3)') nint(1e3*trackrate)
 call compute_wig3j_data_analysis(ell, ellp)
 call cpu_time(tstart)
 if (compute_noise .and. compute_ab_initio) call read_leakage(6,15,ell)
 call cpu_time(tfin)
 print *,'Read time (leakage):',tfin - tstart

 !i = 1
 !call random_seed(size = i)
 !allocate(seeds(i))
 !call random_seed(get=seeds)
 !call random_seed(get=seeds)
 !seed = seeds(size(seeds))

    i = 42
   CALL zigset( i )
 write(lch,'(I3.3)') ell
 write(lch_p,'(I3.3)') ellp

 write(ynum,'(I2.2)') yearnum
 write(ynum2,'(I2.2)') yearnum+nyears-1

 open(102,file='/scratch/shravan/'//instrument//'/tracking'//trackch//'/bcoef_metadata_l_'//lch//'_lp_'//lch_p//&
     '_year_'//ynum//'_'//ynum2,action='write',status='unknown', position='rewind')


 open(112,file='/scratch/shravan/'//instrument//'/tracking'//trackch//'/status_l_'//lch//'_lp_'//lch_p//&
     '_year_'//ynum//'_'//ynum2,action='write',status='unknown', position='rewind')


 if (instrument == 'MDI') then 
   t = 1216 + (yearnum -1)*360
   write(daynum, '(I4.4)') t
   call construct_frequencies ('/scratch/data/'//daynum//'/m10q.36', ell, &
                               ellp,yearnum)
 elseif (instrument == 'HMI') then 
   t = 6328 + (yearnum -1)*360
   write(daynum, '(I4.4)') t
   call construct_frequencies ('/scratch/data/hmi.m10qr.'//daynum//'.36', ell, &
                               ellp,yearnum)

 endif


 allocate(b_coefs(smin:smax,ordmin:ordmax,sigmin:sigmax),tempr(609,305),&
   noise(smin:smax,ordmin:ordmax,sigmin:sigmax),&
   noisevar(smin:smax,ordmin:ordmax,sigmin:sigmax))

 call readfits('/home/shravan/QDP/leaks.fits',tempr,609,305,1)

 do ellt = 0,304
  allocate(leakage(ellt)%ems(-ellt:ellt))
  leakage(ellt)%ems(-ellt:ellt) = tempr(1:2*ellt+1,ellt+1)
 enddo
 deallocate(tempr)


 if (INSTRUMENT == 'MDI') nt = 360 * 24 * 60 * nyears
 if (INSTRUMENT == 'HMI') nt = 360 * 24 * 80 * nyears


  call cpu_time(tstart)

  if (allocated(shtdatap)) deallocate(shtdatap)
  allocate(shtdatap(0:nt-1, 0:ellp))
  shtdatap = shtdat(yearnum, ellp, nyears)

  if (.not. allocated(nu)) then
   allocate(nu(0:nt-1,1))
   if( instrument == 'MDI') nu(0:nt-1,:) = distmat(nt,1)/(nt*60.) * 10**6
   if( instrument == 'HMI') nu(0:nt-1,:) = distmat(nt,1)/(nt*45.) * 10**6
   dnu = nu(2,1) - nu(1,1)
  endif

  if (allocated(shtdata)) deallocate(shtdata)
  allocate(shtdata(0:nt-1, -ell:ell))

  if (ellp == ell) shtdata(:,0:ell) = shtdatap
  if (ellp .ne. ell) shtdata = shtdat(yearnum, ell, nyears)

  do m = -ell,-1
   do offseti = 1, nt-1
    shtdata(offseti,m) = (1-2*modulo(-m,2))*conjg(shtdata(nt-offseti,-m))
   enddo
  enddo
  call cpu_time(tfin)
  print *,'Read time (data):',tfin - tstart

  call cpu_time(tstart)
    !% gammas and nus_res are in muHz
    !% amps are unknown units
   allocate(norms(ordmin:ordmax,-dell:dell))
   norms = 0.d0!1.35e-10/nyears
   do dl = -dell,dell
    eldp = ell + dl
    write(lch_p,'(I3.3)') eldp
    inquire(file='/scratch/shravan/'//instrument//'/norms/'//instrument//&
             '_'//lch_p//'_year_'//ynum,exist=lexist)
    if (lexist) then
     open(1555,file='/scratch/shravan/'//instrument//'/norms/'//instrument//&
             '_'//lch_p//'_year_'//ynum,status='old',action='read')
     do  
      read(1555,*,IOSTAT=ierr) con0, order
      if (ierr .ne.0) exit
      if (order .ge. ordmin .and. order .le. ordmax) &
       norms(order,dl) = con0/nyears
     enddo
     close(1555)
    endif
   enddo
  lch_p = lch

!   open(144,file='/scratch/shravan/'//instrument//&
!     '/norms/'//instrument//'_'//lch//'_year_'//ynum,status='unknown')

   do order = ordmin,ordmax
     do sigoff = sigmin, sigmax
      do s = smin, smax
       allocate(b_coefs(s,order,sigoff)%tee(-s:s),noise(s,order,sigoff)%tee(-s:s),&
              noisevar(s,order,sigoff)%tee(-s:s))
       b_coefs(s,order,sigoff)%tee = 0.d0
       noise(s,order,sigoff)%tee = 0.d0
       noisevar(s,order,sigoff)%tee = 0.d0
      enddo
     enddo
   enddo

   do order = ordmin, ordmax
     if (freqnu(ell)%ords(order) < freqmin) cycle
     minoff = floor((nus(ell)%mn(-ell,order) - 150.d0)/dnu) - 1
     maxoff = floor((nus(ell)%mn(ell,order) + 150.d0)/dnu) + 1

     allocate(prefac(minoff:maxoff,-ell:ell),minoffarr(-ell:ell),tempm(minoff:maxoff),& 
     maxoffarr(-ell:ell),mask(minoff:maxoff,-ell:ell),& !,corr(minoff:maxoff)
     leaks(-delm:delm,-ell:ell),abspow(minoff:maxoff,-ell-dell:ell+dell,-dell:dell), &
     limitpow(minoff:maxoff,-ell-dell:ell+dell,-dell:dell))!,corr(minoff:maxoff))

     abspow = 0.d0
     limitpow = 0.d0
     leaks = 0.d0
     prefac = 0.d0
     nj = 0.d0

     minoffp = minoff
     maxoffp = maxoff

     conr = 0.d0
      do dl = -dell,dell !0,0!
       mix(dl) =0.d0
       nj(dl) = 0.d0
       eldp = ell + dl
       if ((lbound(freqnu(eldp)%ords,1) .gt. order) .or. &
             (ubound(freqnu(eldp)%ords,1) .lt. order)) cycle

       if (eldp > lmax .or. eldp < lmin .or. freqnu(eldp)%ords(order) .lt. freqmin) cycle

       con0 = norms(order,dl)
       if (con0 > 1e-6) & 
        con0 = 0.d0!sum(norms(:,dl))/size(norms(:,dl),1)

!       con0 = 1.d0
       fmode = 278.6/696.d0  * 1e6 * (eldp+0.5) /(4.*pi*2)
       mix(dl) = fmode/freqnu(eldp)%ords(order)**2
       gamnlp = fwhm(eldp)%ords(order)
       nj(dl) = con0 * gamnlp * (amps(eldp)%ords(order) *  & !
       freqnu(eldp)%ords(order)) **2 * (2*pi)**3 * 1e-18

       do emdp =-eldp,eldp
        nulmn = nus(eldp)%mn(emdp,order)
        minoff = max(floor((nulmn - 15*gamnlp)/dnu) -sigmax*3,minoffp)
        maxoff = min(floor((nulmn + 15*gamnlp)/dnu) + sigmax*3,maxoffp)
        limitpow(minoff:maxoff,emdp,dl) = - 1.d0/(nu(minoff:maxoff,1)**2 - &
        (nulmn - cmplx(0.d0,0.5d0)*gamnlp)**2) * 1e12 * contwopi
        abspow(minoff:maxoff,emdp,dl) = abs(limitpow(minoff:maxoff,emdp,dl))**2
       enddo

!       leaks = (rleaks(-delm:delm,:,dl) + mix(dl) * horleaks(-delm:delm,:,dl))**2
!       do m=-ell,ell

!        lolim = floor((nus(ell)%mn(m,order) - 2*gamnlp)/dnu) 
!        hilim = floor((nus(ell)%mn(m,order) + 2*gamnlp)/dnu)

!        do emdp = max(m-delm,-eldp), min(m+delm,eldp)
!         dm = emdp - m
!         conr = sum(abspow(lolim:hilim,emdp,dl))*nj(dl)*leaks(dm,m) + conr
!        enddo
!       enddo

!       if (dl==0) then !print *,nj(dl),freqnu(eldp)%ords(order),gamnlp,'aa',order
!
!         lolim = floor((nus(ell)%mn(0,order) - 2*gamnlp)/dnu) 
!         hilim = floor((nus(ell)%mn(0,order) + 2*gamnlp)/dnu)
!         con = sum(abs(shtdata(lolim:hilim,0)**2.))
!         do m=1,ell
!          
!          lolim = floor((nus(ell)%mn(m,order) - 2*gamnlp)/dnu) 
!          hilim = floor((nus(ell)%mn(m,order) + 2*gamnlp)/dnu)
!
!          con = sum(abs(shtdata(lolim:hilim,m)**2.)) + con
!
!          lolim = floor((nus(ell)%mn(-m,order) - 2*gamnlp)/dnu) 
!          hilim = floor((nus(ell)%mn(-m,order) + 2*gamnlp)/dnu)
!
!          con = sum(abs(shtdata(lolim:hilim,-m)**2.)) + con
!
!          enddo
!        
!       endif

      enddo

!      if (isnan(conr)) con = 0.d0
!      write(144,*) con/conr,order!nus(ell)%mn(m,order),m
!      print *,con,conr,con/conr,order!,nj,freqnu(ell-dell:ell+dell)%ords(order)
!      deallocate(prefac, minoffarr, maxoffarr, limitpow, abspow,leaks,tempm,mask) 
!      cycle

      gamnl = fwhm(ell)%ords(order)

      if (any(en(ell)%ords .eq. order) .and. (freqnu(ellp)%ords(order) .ne. 0)) then

         do t = -smax, smax 

          if (t==0) cycle
          slow = max(smin, abs(t))
          mmin = max(-ell,-ell-t)
          mmax = min(ell,ell-t)

           ! Tracking rate! trackrate
          refind = nint(t*trackrate/dnu) !+ nint(rand())*(1-2*nint(rand()))
          siglow = refind + sigmin !- dsigind
          sighigh = refind + sigmax

          print *,t,order, ell
          write(112,*) t, order, siglow, sighigh, t*a1(ell)%ords(order)
          flush(112)

          do sigoff = siglow,sighigh
            if (sigoff == refind) cycle

            noisehalf = 0.d0
            den(smin:smax) = 0.d0
            sigind = sigoff - refind
            sig = t * trackrate + sigind*dnu 

            do m = mmin, mmax

              mp = m + t
!              if (m==0 .or. (mp==0)) cycle
              nulmn = nus(ell)%mn(m,order)
              nulmnt = nus(ellp)%mn(mp,order)

              lolim_m = nulmn - gamnl
              hilim_m = nulmn + gamnl
              lolim_mp = nulmnt - gamnl - sig
              hilim_mp = nulmnt + gamnl - sig
              
              minoff = floor(min(lolim_m,lolim_mp)/dnu)-1
              maxoff = floor(max(hilim_m,hilim_mp)/dnu)+1 ! mistake, was lolim_m
         
              minoffarr(m) = minoff
              maxoffarr(m) = maxoff

              mask(minoff:maxoff,m) = 0.d0
              where ((((lolim_mp .le. nu(minoff:maxoff,1)) .and. (hilim_mp .ge. nu(minoff:maxoff,1))) &
              .or. ((lolim_m .le. nu(minoff:maxoff,1)) .and. (hilim_m .ge. nu(minoff:maxoff,1))))) mask(minoff:maxoff,m) = 1.d0

              k = 0
              do j=minoff,maxoff
               if (mask(j,m) == 1.d0) then
                k = k+1
                intarr(k) = j
               endif
              enddo

              prefac(intarr(1:k),m) =-2*nu(intarr(1:k),1)*2*pi*1e-6* nj(0)*leakage(ell)%ems(m)*&
              leakage(ellp)%ems(mp)* (limitpow(intarr(1:k)+sigoff,mp,0)*abspow(intarr(1:k),m,0)&
              + conjg(limitpow(intarr(1:k),m,0))*abspow(intarr(1:k)+sigoff,mp,0)) 
             
              conr = sum(abs(prefac(intarr(1:k),m))**2)

              conp =sum(prefac(intarr(1:k),m)*conjg(shtdata(intarr(1:k),m))*shtdata(intarr(1:k)+sigoff,mp)) 
              do s = slow,smax
               b_coefs(s,order, sigind)%tee(t) = b_coefs(s,order,sigind)%tee(t) +  conp*wig(s)%sub3j(t,m)
               den(s) =  conr*wig(s)%sub3j(t,m)**2  + den(s) !
              enddo ! s loop
             enddo ! m loop! omega loop (m')

       !NOISE CALCULATION
           if (compute_noise) then
!             noisevar1 = 0.d0
             do m= mmin,mmax
               minoff = minoffarr(m)


               do dm =0,0!min(mmax-m, delm)!delm)
                emp = m + dm

                maxoff = maxoffarr(emp)

                k = 0
                do j=minoff,maxoff
                 if (mask(j,m)*mask(j,emp) == 1.d0) then
                  k = k+1
                  intarr(k) = j
                 endif
                enddo

                if (k ==0) cycle

                if (compute_ab_initio) then
                 noisehalf = 0.d0
                 noiseother = 0.d0
!                 temphalf = 0.d0
 !                tempother = 0.d0
  
                
                 do dl = -dell,dell
                  eldp = ell + dl
                  leaks = rleaks(-delm:delm,:,dl) + mix(dl) * horleaks(-delm:delm,:,dl)
                  emdpmin = max(-eldp,m-delm,emp-delm)
                  emdpmax = min(eldp,emp+delm,m+delm)

                  do emdp =emdpmin,emdpmax
                   cons = leaks(emdp-m,m) * leaks(emdp-emp,emp)*nj(dl)
                   noisehalf(1:k) = noisehalf(1:k) + cons * abspow(intarr(1:k),emdp,dl)
!                   do reals = 1,nreals
 !                   do indst=1,k
  !                   temphalf(indst,reals) = cons * abspow(intarr(indst),emdp,dl) * &
   !                         (cmplx(rnor(),-rnor()))*cmplx(rnor(),rnor())*0.5 +&
    !                                      temphalf(indst,reals)
     !               enddo
      !             enddo
                  enddo

                  emdpmin = max(m+t-delm,emp+t-delm,-eldp)
                  emdpmax = min(emp+t+delm,m+t+delm,eldp)
             
                  do emdp =emdpmin,emdpmax
                   consp = leaks(emdp-m-t,m+t) * leaks(emdp-emp-t,emp+t)*nj(dl)
                   noiseother(1:k) = noiseother(1:k) + consp * abspow(intarr(1:k)+sigoff,emdp,dl)
!                   do reals = 1,nreals
 !                   do indst=1,k
  !                   tempother(indst,reals) = consp * abspow(intarr(indst)+sigoff,emdp,dl) * &
   !                         (cmplx(rnor(),-rnor()))*cmplx(rnor(),rnor())*0.5 +&
    !                                      tempother(indst,reals)
     ! 
      !              enddo
       !            enddo

                  enddo

                enddo
                conr = (1.d0 + dm/(dm+1e-10))*sum(real(conjg(prefac(intarr(1:k),m))*prefac(intarr(1:k),emp))*&
                               noisehalf(1:k)*noiseother(1:k))
!                do reals =1,nreals
!                  noisevar2(reals) = (1.d0 + dm/(dm+1e-10))*sum(real(conjg(prefac(intarr(1:k),m))*prefac(intarr(1:k),emp)*&
!                      conjg(temphalf(1:k,reals))*tempother(1:k,reals))) - conr
!                enddo

               else
                conr =(1.d0 + dm/(dm+1e-10))*sum(real(conjg(prefac(intarr(1:k),m))*prefac(intarr(1:k),emp)))
               endif

               do s = slow, smax
                cons = wig(s)%sub3j(t,m) * wig(s)%sub3j(t,emp)
                noise(s,order,sigind)%tee(t) = noise(s,order,sigind)%tee(t) + conr * cons 
!                noisevar1(:,s) = noisevar1(:,s) + noisevar2 * cons
               enddo
           
              enddo

             enddo
             
           endif
           den = 1.d0/den
           do s = slow, smax
            b_coefs(s,order,sigind)%tee(t) = b_coefs(s,order,sigind)%tee(t)*den(s)
            if (compute_noise) then
             noise(s,order,sigind)%tee(t) = noise(s,order,sigind)%tee(t) * den(s)**2
!             noisevar(s,order,sigind)%tee(t) = sum(noisevar1(:,s)**2,1)/dble(nreals) * den(s)**4 
            endif
           enddo

       enddo ! sigma loop
      enddo ! t loop

     endif ! to make sure orders represented and the frequency is non-zero

     deallocate(prefac, minoffarr, maxoffarr, limitpow, abspow,leaks,tempm,mask) 
   enddo ! order loop

  call cpu_time(tfin)
  print *,'Computational time:',tfin-tstart
!  call system('rm /scratch/shravan/'//instrument//'/processed/'//ynum//'_'//lch)
!  call exit()
 call writefits_bcoef_same('/scratch/shravan/'//instrument//'/tracking'//trackch//'/', &
      ell, ellp, yearnum, nyears, compute_noise)

end subroutine analyzethis

!================================================================================

function distmat(n,m)
 
 implicit none
 integer m, n, i, j, i2, j2
 real*8 distmat(n,m), sig


 do j=1,m
   j2 = min(j-1,m-j+1)
   do i=1,n
       sig = 1.d0
       i2 = min(i-1,n-i+1)
       if ((i-1) > (n-i+1)) sig = -1.d0
       distmat(i,j) = (i2**2.d0 + j2**2.d0)**0.5d0 * sig
   enddo
 enddo
 
END function distmat

!================================================================================

FUNCTION shtdat(yearnum, ell, nyears)

implicit none

integer *8 fwdplan
integer yearnum, startpoint, nt72, ell, st, twont72
integer indst, indend, nyears, nmax, i, deltast

character*4 daynum,el
character*120 filenam

real*8 norm
real*8, allocatable, dimension(:,:) :: dat

logical lexist
complex*16, allocatable, dimension(:,:) :: inp,shtdat

! 2 years end at 1864 (including it)
!Right after 2 years, i.e. 1936, there is a gap of 180 days and the next dataset
!begins on 2116
!After 2116 there’s another gap of 108 days and the following dataset begins on 2224

!For year 3, start on 2224.
!Ends at 6616 (continuous)


!That’s 61 datasets post vacation and a proper 10 datasets pre-vacation.

!This gives us a total of 70 datasets, i.e. 13 years of data (5 datasets counts
!as a year).

if (instrument == 'MDI') then
  deltast = 0
  startpoint = 1216
  dt = 60.d0
!  domega = 2.d0*pi/(360.d0*24.d0*60.d0*60.d0)
  nt = 360 * 24 * 60 * nyears
  nt72 = 72 * 24 * 60 
  twont72 = 2 * 72 * 24 * 60 
  if (yearnum .gt. 2)  deltast = 288 ! 2224 - 1936

  norm = 1.d0/dble(nt)
  allocate(dat(twont72,ell+1), inp(nt,ell+1), shtdat(nt,ell+1))

  call dfftw_plan_guru_dft(fwdplan,1,nt,1,1,1,ell+1,&
     & nt,nt,inp(1,1),shtdat(1,1),FFTW_BACKWARD,FFTW_ESTIMATE)

  st = (yearnum-1)*360 + startpoint + deltast

  if (ell .ge. 100) write(el,'(I3.3)') ell
  if (ell .lt. 100) write(el,'(I2.2)') ell

  nmax = nyears*5-1

  do i=0,nmax
   write(daynum,'(I4.4)') st
   call readfits('/scratch/data/MDI/mdi.vw_V_sht_gf_72d.'&
       //daynum//'.'//trim(adjustl(el))//'.fits', dat, twont72,ell+1,1)
   indst = i*nt72+1
   indend = (i+1)*nt72
   inp(indst:indend,:) = (dat(1:twont72:2,:) - (0.d0, 1.d0) * dat(2:twont72:2,:))*norm
   st = st + 72
  enddo
!  st = st + 72
!  write(daynum,'(I4.4)') st
!  call readfits('/scratch/data/MDI/mdi.vw_V_sht_gf_72d.'&
!          //daynum//'.'//trim(el)//'.fits', dat, 2*nt72,ell+1,1)
!  inp(nt72+1:2*nt72,:) = (dat(1:2*nt72:2,:) - (0.d0, 1.d0) * dat(2:2*nt72:2,:))*norm
!
!  st = st + 72
!  write(daynum,'(I4.4)') st
!  call readfits('/scratch/data/MDI/mdi.vw_V_sht_gf_72d.'&
!         //daynum//'.'//trim(el)//'.fits', dat, 2*nt72,ell+1,1)
!  inp(2*nt72+1:3*nt72,:) = (dat(1:2*nt72:2,:) - (0.d0, 1.d0) * dat(2:2*nt72:2,:))*norm

!  st = st + 72
!  write(daynum,'(I4.4)') st
!  call readfits('/scratch/data/MDI/mdi.vw_V_sht_gf_72d.'&
!         //daynum//'.'//trim(el)//'.fits', dat, 2*nt72,ell+1,1)
!  inp(3*nt72+1:4*nt72,:) = (dat(1:2*nt72:2,:) - (0.d0, 1.d0) * dat(2:2*nt72:2,:))*norm

!  st = st + 72
!  write(daynum,'(I4.4)') st
!  call readfits('/scratch/data/MDI/mdi.vw_V_sht_gf_72d.'&
!         //daynum//'.'//trim(el)//'.fits', dat, 2*nt72,ell+1,1)
!  inp(4*nt72+1:5*nt72,:) = (dat(1:2*nt72:2,:) - (0.d0, 1.d0) * dat(2:2*nt72:2,:))*norm

  call dfftw_execute_dft(fwdplan,inp,shtdat)

  deallocate(dat,inp)
  call dfftw_destroy_plan(fwdplan)

elseif (instrument == 'HMI') then
  startpoint = 6328
  dt = 45.d0
!  domega = 2.d0*pi/(360.d0*24.d0*60.d0*60.d0)
  nt = 360 * 24 * 80 * nyears
  nt72 = 72 * 24 * 80 
  twont72 = 2 * 72 * 24 * 80 

  norm = 1.d0/dble(nt)
  allocate(dat(twont72,ell+1), inp(nt,ell+1), shtdat(nt,ell+1))

  call dfftw_plan_guru_dft(fwdplan,1,nt,1,1,1,ell+1,&
     & nt,nt,inp(1,1),shtdat(1,1),FFTW_BACKWARD,FFTW_ESTIMATE)

  st = (yearnum-1)*360 + startpoint 

  if (ell .ge. 100) write(el,'(I3.3)') ell
  if (ell .lt. 100) write(el,'(I2.2)') ell

  nmax = nyears*5-1

  do i=0,nmax
   write(daynum,'(I4.4)') st
   call readfits('/scratch/data/HMI/hmi.v_sht_gf_72d.'&
       //daynum//'.'//trim(adjustl(el))//'.fits', dat, twont72,ell+1,1)
   indst = i*nt72+1
   indend = (i+1)*nt72
   inp(indst:indend,:) = (dat(1:twont72:2,:) - (0.d0, 1.d0) * dat(2:twont72:2,:))*norm
   st = st + 72
  enddo

  call dfftw_execute_dft(fwdplan,inp,shtdat)

  deallocate(dat,inp)
  call dfftw_destroy_plan(fwdplan)

endif


END FUNCTION SHTDAT

!================================================================================

        SUBROUTINE writefits_bcoef_same (directory, l, lp, yearnum, nyears,compute_noise)

        implicit none 

        integer blocksize,bitpix,naxes(7),unit1, s, i0, isg
        integer i1, i2, i3, i4, l, lp, yearnum, ind, ns, sig
        integer status1,group,fpixel,flag, nelements, nyears

	character*80 filename
        character*(*) directory
        character*3 lch, lch_p
        character*2 ynum2, ynum
!        character*1 omin, omax

	logical simple,extend, exists, compute_noise

        real*8, allocatable, dimension(:,:,:) :: tempreal
        real*8, allocatable, dimension(:,:,:,:) :: temp


        write(lch,'(I3.3)') l
        write(lch_p,'(I3.3)') lp

!        write(omin,'(I1.1)') ordmin
!        write(omax,'(I1.1)') ordmax
        write(ynum,'(I2.2)') yearnum
        write(ynum2,'(I2.2)') yearnum+nyears-1


         ! WRITING OUT SIGNAL

          filename = directory//'bcoef_l_'//lch//'_lp_'//lch_p//&
          '_year_'//ynum//'_'//ynum2//'.fits'
         
          inquire(file=trim(adjustl(filename)), exist=exists)
          if (exists) call system('rm '//filename)
          print *,'Writing file '//filename

	  status1 = 0
	  call ftgiou(unit1,status1)
	  blocksize=1
!	 dump_array = dble(temp)
	  call ftinit(unit1,trim(adjustl(filename)),blocksize,status1)
	  simple=.true.
	  bitpix=-64
          ns = (smax + 1)**2
	  naxes(1)=ns
	  naxes(2)=ordmax - ordmin + 1
!	  naxes(3)=ordmax - ordmin + 1
	  naxes(3)=sigmax-sigmin+1
	  naxes(4)=2

          allocate(temp(1:ns,ordmin:ordmax,sigmin:sigmax, 2))
          allocate(tempreal(1:ns,ordmin:ordmax,sigmin:sigmax))

          temp = 0.d0
          ind = 0
          do i4=sigmin,sigmax
            do i2 = ordmin, ordmax
             do i1 = smin, smax
              do i0 = -i1, i1
               sig = 1-modulo(i0,2)*2
               ind = i1*(i1+1) + i0 + 1
               temp(ind,i2,i4,1) = &
                 real(b_coefs(i1,i2,i4)%tee(i0))! + &
                     !sig*b_coefs(i1,i2,-i4)%tee(-i0))

               temp(ind,i2,i4,2) = &
                 aimag(b_coefs(i1,i2,i4)%tee(i0))! - &
                     !sig*b_coefs(i1,i2,-i4)%tee(-i0))
              enddo
             enddo
            enddo
          enddo

	  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)!*2
	  extend=.false.
	  group=1
	  fpixel=1

	  call ftphpr(unit1,simple,bitpix,4,naxes,0,1,extend,status1)
	  call ftpprd(unit1,group,fpixel,nelements,temp,status1)
	  call ftclos(unit1, status1)
	  call ftfiou(unit1, status1)

         ! WRITING OUT SIGNAL

          if (compute_noise) then

          filename = directory//'noise_l_'//lch//'_lp_'//lch_p//&
          '_year_'//ynum//'_'//ynum2//'.fits'
         
          inquire(file=trim(adjustl(filename)), exist=exists)
          if (exists) call system('rm '//filename)
          print *,'Writing file '//filename

	  status1 = 0
	  call ftgiou(unit1,status1)
	  blocksize=1
!	 dump_array = dble(temp)
	  call ftinit(unit1,trim(adjustl(filename)),blocksize,status1)
	  simple=.true.
	  bitpix=-64
          ns = (smax + 1)**2
	  naxes(1)=ns
	  naxes(2)=ordmax - ordmin + 1
!	  naxes(3)=ordmax - ordmin + 1
	  naxes(3)=sigmax-sigmin+1

          tempreal = 0.d0
          ind = 0
          do i4=sigmin, sigmax
            do i2 = ordmin, ordmax
             do i1 = smin, smax
              do i0 = -i1, i1
               ind = i1*(i1+1) + i0 + 1
               tempreal(ind,i2,i4) = &
                 real(noise(i1,i2,i4)%tee(i0))! +&
                      !0*noise(i1,i2,-i4)%tee(-i0))
!                 real(noise(i1,i2,i4)%tee(i0))! +&
                      !noise(i1,i2,-i4)%tee(-i0))


              enddo
             enddo
            enddo
          enddo

	  nelements=naxes(1)*naxes(2)*naxes(3)
	  extend=.false.
	  group=1
	  fpixel=1

	  call ftphpr(unit1,simple,bitpix,3,naxes(1:3),0,1,extend,status1)
	  call ftpprd(unit1,group,fpixel,nelements,tempreal,status1)
	  call ftclos(unit1, status1)
	  call ftfiou(unit1, status1)


          if (1 ==2) then 
           filename = directory//'noisevar_l_'//lch//'_lp_'//lch_p//&
           '_year_'//ynum//'_'//ynum2//'.fits'
         
           inquire(file=trim(adjustl(filename)), exist=exists)
           if (exists) call system('rm '//filename)
           print *,'Writing file '//filename

	   status1 = 0
	   call ftgiou(unit1,status1)
	   blocksize=1
!	 dump_array = dble(temp)
	   call ftinit(unit1,trim(adjustl(filename)),blocksize,status1)
	   simple=.true.
	   bitpix=-64
           ns = (smax + 1)**2
	   naxes(1)=ns
	   naxes(2)=ordmax - ordmin + 1
!	  naxes(3)=ordmax - ordmin + 1
	   naxes(3)=sigmax!-sigmin+1

           tempreal = 0.d0
           ind = 0
           do i4=1, sigmax
             do i2 = ordmin, ordmax
              do i1 = smin, smax
               do i0 = -i1, i1
                ind = i1*(i1+1) + i0 + 1
                tempreal(ind,i2,i4) = &
                  real(noisevar(i1,i2,i4)%tee(i0))! +&
                       !(noisevar(i1,i2,-i4)%tee(-i0))

               enddo
              enddo
             enddo
           enddo

	   nelements=naxes(1)*naxes(2)*naxes(3)
	   extend=.false.
	   group=1
	   fpixel=1
 
	   call ftphpr(unit1,simple,bitpix,3,naxes(1:3),0,1,extend,status1)
	   call ftpprd(unit1,group,fpixel,nelements,tempreal,status1)
	   call ftclos(unit1, status1)
	   call ftfiou(unit1, status1)
          
           deallocate(temp,tempreal)
          endif
         endif

         
          do i0 = ordmin, ordmax
            write(102,*) int(en(l)%ords(i0)), freqnu(l)%ords(i0), a1(l)%ords(i0)
          enddo
!           write(102,*) 'Ell: '
!           write(102,*) l
!           write(102,*)
!           write(102,*) 'N: ' 
!           write(102,*) int(en(l)%ords)
!           write(102,*)
!           write(102,*) 'Start year: '
!           write(102,*) yearnum
!           write(102,*)
!           write(102,*) 'End year (including data): '
!           write(102,*) yearnum + nyears - 1
!           write(102,*)
!           write(102,*) 'a_1: '
!           write(102,*) a1(l)%ords(int(en(l)%ords))
!           write(102,*)
!           write(102,*) 'dnu'
!           write(102,*) dnu
           

	 end SUBROUTINE writefits_bcoef_same
!================================================================================


SUBROUTINE compute_wig3j_data_analysis(ell, ellp)

 implicit none
 integer ell, ellp, tmin, tmax, es, em, t, slow, shigh, IER
 real*8 temv(1:(ell+ ellp+ 1)),lmi, lma, con

   temv = 0.d0
   do es = smin, smax
    allocate(wig(es)%sub3j(-es:es,-ell:ell))
    wig(es)%sub3j(-es:es,-ell:ell) = 0.d0
   enddo

   do em = -ell, ell
    
    tmin = max(-smax,-ellp-em) !-em - t <= l'
    tmax = min(smax,-em + ellp) ! -em - t >= -l'
    do t = tmin, tmax

     call drc3jj(dble(ell),dble(ellp),dble(em),dble(-em-t),lmi,lma,temv,lmax+lmax+1,IER) 

      slow = max(int(lmi),smin)
      shigh = min(int(lma),smax)
      con = 1.d0 - 2.d0*modulo(em+t,2)
      do es = slow, shigh
       wig(es)%sub3j(t,em) = temv(es-int(lmi) + 1) * (2.d0*es + 1.d0)**0.5 * con
      enddo
      temv = 0.d0
       
     enddo
    enddo



END SUBROUTINE compute_wig3j_data_analysis


!================================================================================

	 SUBROUTINE readfits(filename,readarr,dim1,dime2,dim3)

	  implicit none
	  integer status,unit,readwrite,blocksize,naxes(3)
	  integer group,firstpix, dim3 ,dime2,dim1
	  integer nelements, hdutype
	  real*8 readarr(dim1,dime2,dim3)
	  real*8 nullval!,temp(dim1,dime2,dim3)
	  logical anynull, lexist
	  character*(*) filename

	      
	  status=0
	  call ftgiou(unit,status)
	  readwrite=0
          inquire(file=filename, exist = lexist)
          if (.not. lexist) then
            print *,filename
            print *,'THIS FILE DOES NOT EXIST'
            stop
          endif
	  print *,'Now reading the file: '//filename


	  call ftopen(unit,filename,readwrite,blocksize,status)

          if (instrument == 'HMI' .and. dim1 > 7e4) &
           call FTMRHD(unit, 1, hdutype, status)

	   naxes(1) = dim1
	   naxes(2) = dime2
	   naxes(3) = dim3
	   nelements=naxes(1)*naxes(2)*naxes(3)
!           print *,nelements
	   group=1
	   firstpix=1
	   nullval=-999

	   call ftgpvd(unit,group,firstpix,nelements,nullval, &
                        readarr,anynull,status)

	   !readarr = temp
	   
!	   print *,minval(readarr),maxval(readarr)
	  call ftclos(unit, status)
	  call ftfiou(unit, status)


	  end SUBROUTINE readfits



!================================================================================

subroutine DRC3JJ (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, &
     IER)
!
!! DRC3JJ evaluates the 3J symbol f(L1) for all allowed values of L1.
!
!***PURPOSE  Evaluate the 3j symbol f(L1) = (  L1   L2 L3)
!                                           (-M2-M3 M2 M3)
!            for all allowed values of L1, the other parameters
!            being held fixed.
!
!***LIBRARY   SLATEC
!***CATEGORY  C19
!***TYPE      DOUBLE PRECISION (RC3JJ-S, DRC3JJ-D)
!***KEYWORDS  3J COEFFICIENTS, 3J SYMBOLS, CLEBSCH-GORDAN COEFFICIENTS,
!             RACAH COEFFICIENTS, VECTOR ADDITION COEFFICIENTS,
!             WIGNER COEFFICIENTS
!***AUTHOR  Gordon, R. G., Harvard University
!           Schulten, K., Max Planck Institute
!***DESCRIPTION
!
! *Usage:
!
!        DOUBLE PRECISION L2, L3, M2, M3, L1MIN, L1MAX, THRCOF(NDIM)
!        INTEGER NDIM, IER
!
!        call DRC3JJ (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, IER)
!
! *Arguments:
!
!     L2 :IN      Parameter in 3j symbol.
!
!     L3 :IN      Parameter in 3j symbol.
!
!     M2 :IN      Parameter in 3j symbol.
!
!     M3 :IN      Parameter in 3j symbol.
!
!     L1MIN :OUT  Smallest allowable L1 in 3j symbol.
!
!     L1MAX :OUT  Largest allowable L1 in 3j symbol.
!
!     THRCOF :OUT Set of 3j coefficients generated by evaluating the
!                 3j symbol for all allowed values of L1.  THRCOF(I)
!                 will contain f(L1MIN+I-1), I=1,2,...,L1MAX+L1MIN+1.
!
!     NDIM :IN    Declared length of THRCOF in calling program.
!
!     IER :OUT    Error flag.
!                 IER=0 No errors.
!                 IER=1 Either L2 < ABS(M2) or L3 < ABS(M3).
!                 IER=2 Either L2+ABS(M2) or L3+ABS(M3) non-integer.
!                 IER=3 L1MAX-L1MIN not an integer.
!                 IER=4 L1MAX less than L1MIN.
!                 IER=5 NDIM less than L1MAX-L1MIN+1.
!
! *Description:
!
!     Although conventionally the parameters of the vector addition
!  coefficients satisfy certain restrictions, such as being integers
!  or integers plus 1/2, the restrictions imposed on input to this
!  subroutine are somewhat weaker. See, for example, Section 27.9 of
!  Abramowitz and Stegun or Appendix C of Volume II of A. Messiah.
!  The restrictions imposed by this subroutine are
!       1. L2  >=  ABS(M2) and L3  >=  ABS(M3);
!       2. L2+ABS(M2) and L3+ABS(M3) must be integers;
!       3. L1MAX-L1MIN must be a non-negative integer, where
!          L1MAX=L2+L3 and L1MIN=MAX(ABS(L2-L3),ABS(M2+M3)).
!  If the conventional restrictions are satisfied, then these
!  restrictions are met.
!
!     The user should be cautious in using input parameters that do
!  not satisfy the conventional restrictions. For example, the
!  the subroutine produces values of
!       f(L1) = ( L1  2.5  5.8)
!               (-0.31.5 -1.2)
!  for L1=3.3,4.3,...,8.3 but none of the symmetry properties of the 3j
!  symbol, set forth on page 1056 of Messiah, is satisfied.
!
!     The subroutine generates f(L1MIN), f(L1MIN+1), ..., f(L1MAX)
!  where L1MIN and L1MAX are defined above. The sequence f(L1) is
!  generated by a three-term recurrence algorithm with scaling to
!  control overflow. Both backward and forward recurrence are used to
!  maintain numerical stability. The two recurrence sequences are
!  matched at an interior point and are normalized from the unitary
!  property of 3j coefficients and Wigner's phase convention.
!
!    The algorithm is suited to applications in which large quantum
!  numbers arise, such as in molecular dynamics.
!
!***REFERENCES  1. Abramowitz, M., and Stegun, I. A., Eds., Handbook
!                  of Mathematical Functions with Formulas, Graphs
!                  and Mathematical Tables, NBS Applied Mathematics
!                  Series 55, June 1964 and subsequent printings.
!               2. Messiah, Albert., Quantum Mechanics, Volume II,
!                  North-Holland Publishing Company, 1963.
!               3. Schulten, Klaus and Gordon, Roy G., Exact recursive
!                  evaluation of 3j and 6j coefficients for quantum-
!                  mechanical coupling of angular momenta, J Math
!                  Phys, v 16, no. 10, October 1975, pp. 1961-1970.
!               4. Schulten, Klaus and Gordon, Roy G., Semiclassical
!                  approximations to 3j  and 6j coefficients for
!                  quantum-mechanical coupling of angular momenta,
!                  J Math Phys, v 16, no. 10, October 1975,
!                  pp. 1971-1988.
!               5. Schulten, Klaus and Gordon, Roy G., Recursive
!                  evaluation of 3j and 6j coefficients, Computer
!                  Phys Comm, v 11, 1976, pp. 269-278.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   880515  SLATEC prologue added by G. C. Nielson, NBS; parameters
!           HUGE and TINY revised to depend on D1MACH.
!   891229  Prologue description rewritten; other prologue sections
!           revised; LMATCH (location of match point for recurrences)
!           removed from argument list; argument IER changed to serve
!           only as an error flag (previously, in cases without error,
!           it returned the number of scalings); number of error codes
!           increased to provide more precise error information;
!           program comments revised; SLATEC error handler calls
!           introduced to enable printing of error messages to meet
!           SLATEC standards. These changes were done by D. W. Lozier,
!           M. A. McClain and J. M. Smith of the National Institute
!           of Standards and Technology, formerly NBS.
!   910415  Mixed type expressions eliminated; variable C1 initialized;
!           description of THRCOF expanded. These changes were done by
!           D. W. Lozier.
!***END PROLOGUE  DRC3JJ
!
  INTEGER NDIM, IER
  DOUBLE PRECISION L2, L3, M2, M3, L1MIN, L1MAX, THRCOF(NDIM)
!
  INTEGER I, INDEX, LSTEP, N, NFIN, NFINP1, NFINP2, NFINP3, NLIM, &
          NSTEP2
  DOUBLE PRECISION A1, A1S, A2, A2S, C1, C1OLD, C2, CNORM, &
                   DENOM, DV, EPS, HUGE, L1, M1, NEWFAC, OLDFAC, &
                   ONE, RATIO, SIGN1, SIGN2, SRHUGE, SRTINY, SUM1, &
                   SUM2, SUMBAC, SUMFOR, SUMUNI, THREE, THRESH, &
                   TINY, TWO, X, X1, X2, X3, Y, Y1, Y2, Y3, ZERO
!
  DATA  ZERO,EPS,ONE,TWO,THREE /0.0D0,0.01D0,1.0D0,2.0D0,3.0D0/
!, D1MACH
!***FIRST EXECUTABLE STATEMENT  DRC3JJ
  IER=0
!  HUGE is the square root of one twentieth of the largest floating
!  point number, approximately.

 !print *,'je',NDIM,L1max-l1min+1
  HUGE = SQRT(D1MACH(2)/20.0D0)
  SRHUGE = SQRT(HUGE)
  TINY = 1.0D0/HUGE
  SRTINY = 1.0D0/SRHUGE
!
!     LMATCH = ZERO
  M1 = - M2 - M3
!
!  Check error conditions 1 and 2.
  if ( (L2-ABS(M2)+EPS < ZERO).OR. &
     (L3-ABS(M3)+EPS < ZERO))THEN
     IER=1
!     call XERMSG('SLATEC','DRC3JJ','L2-ABS(M2) or L3-ABS(M3) '// &
!        'less than zero.',IER,1)
     return
  ELSEIF((MOD(L2+ABS(M2)+EPS,ONE) >= EPS+EPS).OR. &
     (MOD(L3+ABS(M3)+EPS,ONE) >= EPS+EPS))THEN
     IER=2
!     call XERMSG('SLATEC','DRC3JJ','L2+ABS(M2) or L3+ABS(M3) '// &
!        'not integer.',IER,1)
     return
  end if
!
!
!
!  Limits for L1
!
!  print *,l2,l3,m1,l1min,l1max,max(abs(l2-l3),abs(m1))
!  stop
  L1MIN = MAX(ABS(L2-L3),ABS(M1))!, dble(smin))
  L1MAX = (L2 + L3)!, dble(smax))

!  Check error condition 3.
  if ( MOD(L1MAX-L1MIN+EPS,ONE) >= EPS+EPS)THEN
     IER=3
 !    call XERMSG('SLATEC','DRC3JJ','L1MAX-L1MIN not integer.',IER,1)
     return
  end if
  if ( L1MIN < L1MAX-EPS)   go to 20
  if ( L1MIN < L1MAX+EPS)   go to 10
!
!  Check error condition 4.
  IER=4
 ! call XERMSG('SLATEC','DRC3JJ','L1MIN greater than L1MAX.',IER,1)
  return
!
!  This is reached in case that L1 can take only one value,
!  i.e. L1MIN = L1MAX
!
   10 CONTINUE

!     LSCALE = 0
  THRCOF(1) = (-ONE) ** INT(ABS(L2+M2-L3+M3)+EPS) / &
   SQRT(L1MIN + L2 + L3 + ONE)
  return
!
!  This is reached in case that L1 takes more than one value,
!  i.e. L1MIN < L1MAX.
!
   20 CONTINUE
!     LSCALE = 0
  NFIN = INT(L1MAX-L1MIN+ONE+EPS)
  if ( NDIM-NFIN)  21, 23, 23
!
!  Check error condition 5.
   21 IER = 5
 ! call XERMSG('SLATEC','DRC3JJ','Dimension of result array for '// &
 !             '3j coefficients too small.',IER,1)
  return
!
!
!  Starting forward recursion from L1MIN taking NSTEP1 steps
!
   23 L1 = L1MIN
  NEWFAC = 0.0D0
  C1 = 0.0D0
  THRCOF(1) = SRTINY
  SUM1 = (L1+L1+ONE) * TINY
!
!
  LSTEP = 1
   30 LSTEP = LSTEP + 1
  L1 = L1 + ONE
!
!
  OLDFAC = NEWFAC
  A1 = (L1+L2+L3+ONE) * (L1-L2+L3) * (L1+L2-L3) * (-L1+L2+L3+ONE)
  A2 = (L1+M1) * (L1-M1)
  NEWFAC = SQRT(A1*A2)
  if ( L1 < ONE+EPS)   go to 40
!
!
  DV = - L2*(L2+ONE) * M1 + L3*(L3+ONE) * M1 + L1*(L1-ONE) * (M3-M2)
  DENOM = (L1-ONE) * NEWFAC
!
  if ( LSTEP-2)  32, 32, 31
!
   31 C1OLD = ABS(C1)
   32 C1 = - (L1+L1-ONE) * DV / DENOM
  go to 50
!
!  If L1 = 1, (L1-1) has to be factored out of DV, hence
!
   40 C1 = - (L1+L1-ONE) * L1 * (M3-M2) / NEWFAC
!
   50 if ( LSTEP > 2)   go to 60
!
!
!  If L1 = L1MIN + 1, the third term in the recursion equation vanishes,
!  hence
  X = SRTINY * C1
  THRCOF(2) = X
  SUM1 = SUM1 + TINY * (L1+L1+ONE) * C1*C1
  if ( LSTEP == NFIN)   go to 220
  go to 30
!
!
   60 C2 = - L1 * OLDFAC / DENOM
!
!  Recursion to the next 3j coefficient X
!
  X = C1 * THRCOF(LSTEP-1) + C2 * THRCOF(LSTEP-2)
  THRCOF(LSTEP) = X
  SUMFOR = SUM1
  SUM1 = SUM1 + (L1+L1+ONE) * X*X
  if ( LSTEP == NFIN)   go to 100
!
!  See if last unnormalized 3j coefficient exceeds SRHUGE
!
  if ( ABS(X) < SRHUGE)   go to 80
!
!  This is reached if last 3j coefficient larger than SRHUGE,
!  so that the recursion series THRCOF(1), ... , THRCOF(LSTEP)
!  has to be rescaled to prevent overflow
!
!     LSCALE = LSCALE + 1
  DO 70 I=1,LSTEP
  if ( ABS(THRCOF(I)) < SRTINY)   THRCOF(I) = ZERO
   70 THRCOF(I) = THRCOF(I) / SRHUGE
  SUM1 = SUM1 / HUGE
  SUMFOR = SUMFOR / HUGE
  X = X / SRHUGE
!
!  As long as ABS(C1) is decreasing, the recursion proceeds towards
!  increasing 3j values and, hence, is numerically stable.  Once
!  an increase of ABS(C1) is detected, the recursion direction is
!  reversed.
!
   80 if ( C1OLD-ABS(C1))   100, 100, 30
!
!
!  Keep three 3j coefficients around LMATCH for comparison with
!  backward recursion.
!
  100 CONTINUE
!     LMATCH = L1 - 1
  X1 = X
  X2 = THRCOF(LSTEP-1)
  X3 = THRCOF(LSTEP-2)
  NSTEP2 = NFIN - LSTEP + 3
!
!
!
!
!  Starting backward recursion from L1MAX taking NSTEP2 steps, so
!  that forward and backward recursion overlap at three points
!  L1 = LMATCH+1, LMATCH, LMATCH-1.
!
  NFINP1 = NFIN + 1
  NFINP2 = NFIN + 2
  NFINP3 = NFIN + 3
  L1 = L1MAX
  THRCOF(NFIN) = SRTINY
  SUM2 = TINY * (L1+L1+ONE)
!
  L1 = L1 + TWO
  LSTEP = 1
  110 LSTEP = LSTEP + 1
  L1 = L1 - ONE
!
  OLDFAC = NEWFAC
  A1S = (L1+L2+L3)*(L1-L2+L3-ONE)*(L1+L2-L3-ONE)*(-L1+L2+L3+TWO)
  A2S = (L1+M1-ONE) * (L1-M1-ONE)
  NEWFAC = SQRT(A1S*A2S)
!
  DV = - L2*(L2+ONE) * M1 + L3*(L3+ONE) * M1 + L1*(L1-ONE) * (M3-M2)
!
  DENOM = L1 * NEWFAC
  C1 = - (L1+L1-ONE) * DV / DENOM
  if ( LSTEP > 2)   go to 120
!
!  If L1 = L1MAX + 1, the third term in the recursion formula vanishes
!
  Y = SRTINY * C1
  THRCOF(NFIN-1) = Y
  SUMBAC = SUM2
  SUM2 = SUM2 + TINY * (L1+L1-THREE) * C1*C1
!
  go to 110
!
!
  120 C2 = - (L1 - ONE) * OLDFAC / DENOM
!
!  Recursion to the next 3j coefficient Y
!
  Y = C1 * THRCOF(NFINP2-LSTEP) + C2 * THRCOF(NFINP3-LSTEP)
!
  if ( LSTEP == NSTEP2)   go to 200
!
  THRCOF(NFINP1-LSTEP) = Y
  SUMBAC = SUM2
  SUM2 = SUM2 + (L1+L1-THREE) * Y*Y
!
!  See if last unnormalized 3j coefficient exceeds SRHUGE
!
  if ( ABS(Y) < SRHUGE)   go to 110
!
!  This is reached if last 3j coefficient larger than SRHUGE,
!  so that the recursion series THRCOF(NFIN), ... ,THRCOF(NFIN-LSTEP+1)
!  has to be rescaled to prevent overflow
!
!     LSCALE = LSCALE + 1
  DO 130 I=1,LSTEP
  INDEX = NFIN - I + 1
  if ( ABS(THRCOF(INDEX)) < SRTINY)   THRCOF(INDEX) = ZERO
  130 THRCOF(INDEX) = THRCOF(INDEX) / SRHUGE
  SUM2 = SUM2 / HUGE
  SUMBAC = SUMBAC / HUGE
!
!
  go to 110
!
!
!  The forward recursion 3j coefficients X1, X2, X3 are to be matched
!  with the corresponding backward recursion values Y1, Y2, Y3.
!
  200 Y3 = Y
  Y2 = THRCOF(NFINP2-LSTEP)
  Y1 = THRCOF(NFINP3-LSTEP)
!
!
!  Determine now RATIO such that YI = RATIO * XI  (I=1,2,3) holds
!  with minimal error.
!
  RATIO = ( X1*Y1 + X2*Y2 + X3*Y3 ) / ( X1*X1 + X2*X2 + X3*X3 )
  NLIM = NFIN - NSTEP2 + 1
!
  if ( ABS(RATIO) < ONE)   go to 211
!
  DO 210 N=1,NLIM
  210 THRCOF(N) = RATIO * THRCOF(N)
  SUMUNI = RATIO * RATIO * SUMFOR + SUMBAC
  go to 230
!
  211 NLIM = NLIM + 1
  RATIO = ONE / RATIO
  DO 212 N=NLIM,NFIN
  212 THRCOF(N) = RATIO * THRCOF(N)
  SUMUNI = SUMFOR + RATIO*RATIO*SUMBAC
  go to 230
!
  220 SUMUNI = SUM1
!
!
!  Normalize 3j coefficients
!
  230 CNORM = ONE / SQRT(SUMUNI)
!
!  Sign convention for last 3j coefficient determines overall phase
!
  SIGN1 = SIGN(ONE,THRCOF(NFIN))
  SIGN2 = (-ONE) ** INT(ABS(L2+M2-L3+M3)+EPS)
  if ( SIGN1*SIGN2) 235,235,236
  235 CNORM = - CNORM
!
  236 if ( ABS(CNORM) < ONE)   go to 250
!
  DO 240 N=1,NFIN
  240 THRCOF(N) = CNORM * THRCOF(N)
  return
!
  250 THRESH = TINY / ABS(CNORM)
  DO 251 N=1,NFIN
  if ( ABS(THRCOF(N)) < THRESH)   THRCOF(N) = ZERO
  251 THRCOF(N) = CNORM * THRCOF(N)
!
  return
end subroutine DRC3JJ

!================================================================================

FUNCTION D1MACH (I)
!
!! D1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        D = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of D above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH( 3) = B**(-T), the smallest relative spacing.
!   D1MACH( 4) = B**(1-T), the largest relative spacing.
!   D1MACH( 5) = LOG10(B)
!
!   Assume double precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(14) = T, the number of base-B digits.
!   I1MACH(15) = EMIN, the smallest exponent E.
!   I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   900911  Added SUN 386i constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  D1MACH
!
  double precision d1mach
  INTEGER I
  INTEGER SMALL(4)
  INTEGER LARGE(4)
  INTEGER RIGHT(4)
  INTEGER DIVER(4)
  INTEGER LOG10(4)
!
  DOUBLE PRECISION DMACH(5)
  SAVE DMACH
!
  EQUIVALENCE (DMACH(1),SMALL(1))
  EQUIVALENCE (DMACH(2),LARGE(1))
  EQUIVALENCE (DMACH(3),RIGHT(1))
  EQUIVALENCE (DMACH(4),DIVER(1))
  EQUIVALENCE (DMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
!     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
!     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
!     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA SMALL(1) / ZC00800000 /
!     DATA SMALL(2) / Z000000000 /
!     DATA LARGE(1) / ZDFFFFFFFF /
!     DATA LARGE(2) / ZFFFFFFFFF /
!     DATA RIGHT(1) / ZCC5800000 /
!     DATA RIGHT(2) / Z000000000 /
!     DATA DIVER(1) / ZCC6800000 /
!     DATA DIVER(2) / Z000000000 /
!     DATA LOG10(1) / ZD00E730E7 /
!     DATA LOG10(2) / ZC77800DC0 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O0000000000000000 /
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O0007777777777777 /
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O7770000000000000 /
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O7777777777777777 /
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA SMALL(1) / Z"3001800000000000" /
!     DATA SMALL(2) / Z"3001000000000000" /
!     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
!     DATA LARGE(2) / Z"4FFE000000000000" /
!     DATA RIGHT(1) / Z"3FD2800000000000" /
!     DATA RIGHT(2) / Z"3FD2000000000000" /
!     DATA DIVER(1) / Z"3FD3800000000000" /
!     DATA DIVER(2) / Z"3FD3000000000000" /
!     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
!     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA SMALL(1) / 00564000000000000000B /
!     DATA SMALL(2) / 00000000000000000000B /
!     DATA LARGE(1) / 37757777777777777777B /
!     DATA LARGE(2) / 37157777777777777777B /
!     DATA RIGHT(1) / 15624000000000000000B /
!     DATA RIGHT(2) / 00000000000000000000B /
!     DATA DIVER(1) / 15634000000000000000B /
!     DATA DIVER(2) / 00000000000000000000B /
!     DATA LOG10(1) / 17164642023241175717B /
!     DATA LOG10(2) / 16367571421742254654B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn OR -pd8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CC0000000000000' /
!     DATA DMACH(4) / Z'3CD0000000000000' /
!     DATA DMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
!     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
!     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA SMALL(1) / 201354000000000000000B /
!     DATA SMALL(2) / 000000000000000000000B /
!     DATA LARGE(1) / 577767777777777777777B /
!     DATA LARGE(2) / 000007777777777777774B /
!     DATA RIGHT(1) / 376434000000000000000B /
!     DATA RIGHT(2) / 000000000000000000000B /
!     DATA DIVER(1) / 376444000000000000000B /
!     DATA DIVER(2) / 000000000000000000000B /
!     DATA LOG10(1) / 377774642023241175717B /
!     DATA LOG10(2) / 000007571421742254654B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC DMACH(5)
!
!     DATA SMALL /    20K, 3*0 /
!     DATA LARGE / 77777K, 3*177777K /
!     DATA RIGHT / 31420K, 3*0 /
!     DATA DIVER / 32020K, 3*0 /
!     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA DMACH(1) / '0000000000000010'X /
!     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
!     DATA DMACH(3) / '0000000000003CC0'X /
!     DATA DMACH(4) / '0000000000003CD0'X /
!     DATA DMACH(5) / '79FF509F44133FF3'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FORMAT
!
      DATA DMACH(1) / Z'0010000000000000' /
      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      DATA DMACH(3) / Z'3CA0000000000000' /
      DATA DMACH(4) / Z'3CB0000000000000' /
      DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
!     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
!     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
!     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
!     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1), SMALL(2) /        128,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!     DATA DIVER(1), DIVER(2) /       9472,           0 /
!     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
!
!     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1), SMALL(2) /         16,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!     DATA DIVER(1), DIVER(2) /      15568,           0 /
!     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
!
!     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
!
!     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
!     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
!     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
!     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
!     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
!     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
!     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
!     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
!     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) /  40000B,       0 /
!     DATA SMALL(3), SMALL(4) /       0,       1 /
!     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
!     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
!     DATA RIGHT(3), RIGHT(4) /       0,    225B /
!     DATA DIVER(1), DIVER(2) /  40000B,       0 /
!     DATA DIVER(3), DIVER(4) /       0,    227B /
!     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
!     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
!     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
!     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
!     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
!     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
!     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
!     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
!
!     DATA SMALL(1) / 2.23D-308  /
!     DATA LARGE(1) / 1.79D+308  /
!     DATA RIGHT(1) / 1.11D-16   /
!     DATA DIVER(1) / 2.22D-16   /
!     DATA LOG10(1) / 0.301029995663981195D0 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
!     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
!     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /    8388608,           0 /
!     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
!     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
!     DATA DIVER(1), DIVER(2) /  620756992,           0 /
!     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
!
!     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
!     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
!     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
!     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
!     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /    128,      0 /
!     DATA SMALL(3), SMALL(4) /      0,      0 /
!     DATA LARGE(1), LARGE(2) /  32767,     -1 /
!     DATA LARGE(3), LARGE(4) /     -1,     -1 /
!     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
!     DATA RIGHT(3), RIGHT(4) /      0,      0 /
!     DATA DIVER(1), DIVER(2) /   9472,      0 /
!     DATA DIVER(3), DIVER(4) /      0,      0 /
!     DATA LOG10(1), LOG10(2) /  16282,   8346 /
!     DATA LOG10(3), LOG10(4) / -31493, -12296 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA SMALL(3), SMALL(4) / O000000, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA LARGE(3), LARGE(4) / O177777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
!     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
!     DATA DIVER(1), DIVER(2) / O022400, O000000 /
!     DATA DIVER(3), DIVER(4) / O000000, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020232 /
!     DATA LOG10(3), LOG10(4) / O102373, O147770 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
!     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
!     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
!
!     MACHINE CONSTANTS FOR THE SUN 386i
!
!     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
!     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
!     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
!     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
!
  if ( I < 1 .OR. I > 5 ) then
!    call XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

  D1MACH = DMACH(I)

  return
end function D1MACH

!================================================================================================
subroutine read_leakage(dl,dm, ell)
  implicit none
  integer dm, dl, dl_mat,dm_mat, ind, ind0, i, j, l, lp, m, mp, ell
  real*8 L_rr,L_ii
  ! integer dm,dl
  !integer, allocatable, dimension(:) :: lp_mat,mp_mat
  !integer numlines,line_no_match,i,unit,io
  real*8, dimension(:,:,:),allocatable :: cr_mat,ci_mat, hr_mat, hi_mat

  dl_mat=6
  dm_mat=15

!  if (abs(lp -l) .le. dl_mat .and. (abs(mp-m) .le. dm_mat)) then
     

  ! dl = l - lp
  ! dm = m - mp

!   unit=100

!   open(unit,FILE="leakvw0/default/list320")
!   numlines = 0
!   do
!     read(unit,*,iostat=io)
!     if (io/=0) exit
!     numlines = numlines+1
!   enddo

!   numlines = 51681
!   allocate(lp_mat(numlines),mp_mat(numlines))

 !  rewind(unit)

!   do i=1,numlines
!     read(unit,*) lp_mat(i),mp_mat(i)
!     if ((lp_mat(i)==lp) .and. (mp_mat(i)==mp)) then
!       line_no_match=i
!       exit
!     endif
!   enddo

!   close(unit)

   ind = 51681 !(321 * 322 / 2)
   allocate(cr_mat(2*dm_mat+1,2*dl_mat+1,ind),ci_mat(2*dm_mat+1,2*dl_mat+1,ind), &
         hr_mat(2*dm_mat+1,2*dl_mat+1,ind),hi_mat(2*dm_mat+1,2*dl_mat+1,ind))

   call readfits('/home/shravan/QDP/leakvw0/default/vradsum/leakrlist.vradsum.fits',cr_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)
   call readfits('/home/shravan/QDP/leakvw0/default/vradsum/leakilist.vradsum.fits',ci_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)

   call readfits('/home/shravan/QDP/leakvw0/default/vhorsum/leakrlist.vhorsum.fits',hr_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)
   call readfits('/home/shravan/QDP/leakvw0/default/vhorsum/leakilist.vhorsum.fits',hi_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)

   
   allocate(rleaks(-dm:dm,-ell:ell,-dl:dl), horleaks(-dm:dm,-ell:ell,-dl:dl))
   rleaks = 0.d0
   horleaks = 0.d0
   do lp = ell,ell!lmin, lmax
    ind0 = lp*(lp+1)/2
!    do j = -dl, dl
!     do i = -dm, dm
!      allocate(rleaks(i,j,lp)%em(-lp:lp),horleaks(i,j,lp)%em(-lp:lp))
!     enddo
!    enddo

    do mp = 0, lp
     ind = ind0 + mp + 1
     do j = -dl, dl
      l = lp + j

      do i = -dm, dm

       m = mp + i
       if (abs(m) > l) cycle

       L_ii = 0.5*sign(1,m)*ci_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)
       L_rr = 0.5*cr_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)

       !if(abs(L_rr*L_ii) > 1e-2) cycle

       rleaks(i,mp,j) = L_rr + L_ii  
       rleaks(i,-mp,j) = L_rr + sign(1,mp)* L_ii

       L_ii = 0.5*sign(1,m)*hi_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)
       L_rr = 0.5*hr_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)

       horleaks(i,mp,j) = L_rr + L_ii  
       horleaks(i,-mp,j) = L_rr + sign(1,mp)* L_ii

!       if(abs(L_rr*L_ii) > 1e-2)   print *,rleaks(i,mp,j),horleaks(i,mp,j),m,mp,ell,lp,dm
       !   rleaks(i,mp,j) = 0.d0 !
       !   horleaks(i,mp,j) = 0.d0 ! print *,rleaks(i,mp,j),horleaks(i,mp,j),m,mp,ell,lp,dm
       !endif
      enddo
     enddo
    enddo
   enddo

   deallocate(cr_mat,ci_mat,hr_mat,hi_mat)
!stop
   !leak = (cr+ci)/2.d0

!  else

   !leak = 0.d0
!   print *,'BEYOND ACCEPTABLE RANGE!'

!  endif

 

end subroutine read_leakage



!================================================================================

!C *************************************************************************
      subroutine readheader

!C  Print out all the header keywords in all extensions of a FITS file

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      character filename*80,record*80

!C  The STATUS parameter must always be initialized.
      status=0

!C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

!C     name of FITS file 
      filename='/scratch/data/HMI/hmi.v_sht_gf_72d.6472.120.fits'

!C     open the FITS file, with read-only access.  The returned BLOCKSIZE
!C     parameter is obsolete and should be ignored. 
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

      j = 0
100   continue
      j = j + 1

      print *,'Header listing for HDU', j

!C  The FTGHSP subroutine returns the number of existing keywords in the
!C  current header data unit (CHDU), not counting the required END keyword,
      call ftghsp(unit,nkeys,nspace,status)

!C  Read each 80-character keyword record, and print it out.
      do i = 1, nkeys
          call ftgrec(unit,i,record,status)
          print *,record
      end do

!C  Print out an END record, and a blank line to mark the end of the header.
      if (status .eq. 0)then
          print *,'END'
          print *,' '
      end if

!C  Try moving to the next extension in the FITS file, if it exists.
!C  The FTMRHD subroutine attempts to move to the next HDU, as specified by
!C  the second parameter.   This subroutine moves by a relative number of
!C  HDUs from the current HDU.  The related FTMAHD routine may be used to
!C  move to an absolute HDU number in the FITS file.  If the end-of-file is
!C  encountered when trying to move to the specified extension, then a
!C  status = 107 is returned.
      call ftmrhd(unit,1,hdutype,status)

      if (status .eq. 0)then
!C         success, so jump back and print out keywords in this extension
          go to 100

      else if (status .eq. 107)then
!C         hit end of file, so quit
          status=0
      end if

!C  The FITS file must always be closed before exiting the program. 
!C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

!C  Check for any error, and if so print out error messages.
!C  The PRINTERROR subroutine is listed near the end of this file.
      !if (status .gt. 0)call printerror(status)
      end subroutine readheader
!C *************************************************************************
END MODULE DATA_ANALYSIS