MODULE KERNELS

use frequencies
implicit none
!include 'params.i'
include 'fftw3.f'

integer,allocatable, dimension(:,:) :: mapfwd
real*8 freqnorm, rhonorm, speednorm !pi, rsun, msun, Gcon, 
!parameter (pi = acos(-1.0d0), Gcon = 6.678e-8)
real*8, allocatable, dimension(:,:) :: eigU, eigV, deigU, deigV
real*8, allocatable, dimension(:) :: omegaref, r, rho, &
                    stretch, unstretch, dr, drhodr, c2
complex*16, allocatable, dimension(:,:) :: w_st
integer orders(lmin:lmax,2)
   
type :: wigcell
  real*8, allocatable, dimension(:,:) :: sub3j
end type  
 
type(wigcell) :: wigner(smin:smax,lmin:lmax,lmin:lmax), wig(smin:smax)
type(wigcell), allocatable, dimension(:,:,:) :: gammafun, gammasptp

type :: ems
  real*8, allocatable, dimension(:) :: em, sub3j, en
end type  


type(ems),allocatable, dimension(:) :: gammast
type(ems),allocatable, dimension(:,:) :: leakmat
type(ems), allocatable, dimension(:,:,:) :: rleaks, horleaks
type(ems), allocatable, dimension(:) :: freqsbase, leakage
real*8, allocatable, dimension(:,:,:,:,:) :: Qmat

logical existence(lmin:lmax,0:30)
 logical compute_wigner

Contains

!------------------

subroutine basic_setup!(compute_wigner)

 implicit none
! include 'params.i'

 integer k, em, t, ell, ier, indi, es, ellp
 integer tmin, tmax, slow, shigh, io, n
 real*4 arr(100)
 real*8 temv(1:(lmax+lmax+1)), lmi, lma
 real*8, dimension(:), allocatable :: tempa
 real*8, dimension(:,:), allocatable :: tempr

   !sigmin = int(7.2*8.64*0.25)
   !sigmax = int(7.2*8.64*6)
   ordmax = 30
   ordmin = 0
   open(99, file='/scratch/jb6888/QDP/egvt.sfopal5h5', form='unformatted', status='old')
   read(99)
   read(99) arr
   close(99)
   nr = floor(arr(21)) -1
!   ntot = (lmax-lmin+1) * (ordmax - ordmin+1)
   ns = smax - smin + 1


!   open(44,file='freqs_base',position='rewind',action='read')
!   allocate(freqsbase(lmin:lmax))

!   do k=lmin,lmax
!    allocate(freqsbase(k)%en(0:24))
!   enddo
 
!   ntot = 0
!   do
!     read(44,*,iostat=io) ell, n, lmi, lma
!     if (io < 0) exit
!     if (ell .ge. lmin .and. ell .le. lmax) then
!      freqsbase(int(ell))%en(int(n)) = dble(lmi)
!!      existence(ell,n) = .true.
!      ntot = ntot + 1
!     endif
!   enddo
!   close(44)
  !if (instrument == 'MDI')
  if (instrument == 'MDI') call construct_frequencies(trim(freqmdidir)//'/mdi.1216.36',lmin,lmax,1)
  if (instrument == 'HMI')call construct_frequencies(trim(freqhmidir)//'/hmi.6328.36',lmin,lmax,1)

  ntot = 0
  existence(:,:) =.false.
  orders(:,1) = 45
  orders(:,2) = 0

  do ell = lmin,lmax
   do n = lbound(en(ell)%ords,1),ubound(en(ell)%ords,1)
    if (freqnu(ell)%ords(n) .ne. 0.d0 .and. en(ell)%ords(n) .ne. -1.d0) then
      existence(ell,n) = .true.
      if (n .le. orders(int(ell),1)) orders(int(ell),1) = n
      if (n .ge. orders(int(ell),2)) orders(int(ell),2) = n
      ntot = ntot + 1
    endif
   enddo
  enddo

   allocate(r(nr), rho(nr), stretch(nr), unstretch(nr), &
   eigU(nr,ntot),deigU(nr,ntot),eigV(nr,ntot), omegaref(ntot), &
   deigV(nr,ntot),mapfwd(lmin:lmax,0:30), dr(nr), drhodr(nr),&
   tempa(nr), c2(nr))


   call binaryreader ()
 
   call dbyd1(unstretch,r,nr,1,1)

   do k=1,nr
    stretch(k) = 1.0/unstretch(k)
   enddo

   tempa(1:nr) = r(1:nr)*log(rho(1:nr))
   call dbyd1(deigU, eigU, nr, ntot, 1)
   call dbyd1(deigV, eigV, nr, ntot, 1)
   call dbyd1(drhodr, tempa, nr, 1, 1)

   deallocate(tempa)
   do k=1,ntot
    deigU(:,k) = deigU(:,k) * stretch
    deigV(:,k) = deigV(:,k) * stretch
   enddo
   drhodr = drhodr * stretch

!   compute_wigner = .true.
   if (compute_wigner) then
   temv = 0.d0
   do ellp = lmin, lmax
   ! do ell = lmin, lmax
    ell = ellp

     do es = smin, smax
      allocate(wigner(es,ell,ellp)%sub3j(-es:es,-ell:ell))
      wigner(es,ell,ellp)%sub3j(-es:es,-ell:ell) = 0.d0
     enddo

     do em = -ell, ell
    
      tmin = max(-smax,-ellp-em) !-em - t <= l'
      tmax = min(smax,-em + ellp) ! -em - t >= -l'
      !print *,ell,ellp,em,tmin,tmax
      do t = tmin, tmax!-smax, smax 

       call drc3jj(dble(ell),dble(ellp),dble(em),dble(-em-t),lmi,lma,temv,lmax+lmax+1,IER) 

!       do indi = 1,int(lma-lmi+1)
!        if (abs(temv(indi)) .ge. 1e4) temv = 0.d0
!        if (abs(temv(indi)) .le. 1e-19) temv = 0.d0
!       enddo
       !if (ier .ne. 0) print *,-em - t,lmi,lma,ier
!       if (IER == 0) then
        slow = max(int(lmi),smin)
        shigh = min(int(lma),smax)

        do es = slow, shigh
         !print *,lmi, t, tmin, tmax,temv(es-int(lmi)+1), em,es,ell
         wigner(es,ell,ellp)%sub3j(t,em) = temv(es-int(lmi) + 1)
!         if (t==0 .and. es==21) print *,slow,shigh,ell,ellp,em,temv(es-int(lmi)+1),(lmi),lma, IER

!	if (ellp == 100) print *,es,ell,ellp,slow,shigh,wigner(es,ell,ellp)%sub3j(es-5,ell-5)
        enddo
!        print *,wigner(slow,ell,ellp)%sub3j(slow,em),slow,em
        temv = 0.d0
       
!       endif

       !if (em == -10 .and. int(lmi) .le. 60) print *,int(lmi),t,wigner(ell)%sub3j(int(lmi),em,t),temv(1)
!       if (em == -20 .and. t == -10) then
!           do k=1,smax
!            print *,wigner(ell)%sub3j(k,em,t),k,em,t
!           enddo
!        endif
         ! print *,lma,lmi
    !     do k=1,int(lma-lmi+1)!,int(lma)
    !      if (k+int(lmi)-1 .le. smax) &
    !      print *, temv(k),k+int(lmi)-1!lmi,lma
    !     enddo
       !endif

       !print *,em,t,maxval(abs(temv(1:int(lma-lmi+2))))
!     enddo
    enddo
   enddo
 enddo
 endif
        !print *,wigner(4,120,120)%sub3j(-5,116)!,slow,em

! open(443,file='wigner_comparison',action='write',status='replace')
! ell = 120
 !es = 20
! em = 30
 
! do es = 4,35 !em = -ell, ell
!  write(443,*) es,wigner(es,ell,ell)%sub3j(es-3:es,em)
! enddo
! close(443)
 !print *,size(wigner(es-1,ell,ell)%sub3j),(2*es-1)*(241)
! stop

 allocate(tempr(609,305),leakage(0:304))

 call readfits('/scratch/jb6888/QDP/leaks.fits',tempr,609,305,1)
 do ell = 0,304
  allocate(leakage(ell)%em(-ell:ell))
  leakage(ell)%em(-ell:ell) = tempr(1:2*ell+1,ell+1)
 enddo
 deallocate(tempr)

end subroutine basic_setup


!------------------

subroutine compute_kernels_all(s, dell, asymptotics)  !inversion !(polR, polH, polT, s, leng)

implicit none
integer*8 planfwd, planinv
integer s, ind, indp, leng, ord, k, ell, ellp, ordp,jj, dell,storderp
!parameter(leng = (lmax-lmin+1)**2*(ordmax-ordmin+1)**2)
complex*16, dimension(:,:), allocatable :: polR, polH, polT, polc
real*8 conl, conlp, overallcon, l, m, n, lp, mp, np
real*8 conomeg, con2, tempa(nr), outp(nr), dr(nr), freqdiff
complex*16 coeffs(1:floor(nr/2.d0)+1), filter(1:floor(nr/2.d0)+1)
logical asymptotics
character*3 llow, lhigh
character*1 dellst
 
dr(1:nr-1) = r(2:nr)-r(1:nr-1)
dr(nr) = dr(nr-1)



 call dfftw_plan_dft_r2c_1d(planfwd,nr,tempa(1), coeffs(1), FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_1d(planinv,nr,coeffs(1), tempa(1), FFTW_ESTIMATE)

coeffs = 0.d0
do k=1,floor(nr/2.d0)+1
 filter(k) = cmplx(1.d0/(1.d0 + exp(k - 1 - nr/30.d0)), 0.d0)/dble(nr)
enddo

leng=0
 do ell = lmin, lmax
  ellp = ell + dell
  do ord = orders(ell,1), orders(ell,2)
   storderp = lbound(en(ellp)%ords,1)
   if (ell == ellp) storderp = ord
   do ordp = storderp, orders(ellp,2)
    if (existence(ell, ord) .and. existence(ellp,ordp)) then
     freqdiff = abs(freqnu(ell)%ords(ord) - freqnu(ellp)%ords(ordp))
     if (freqdiff .gt. sigmax .or. freqdiff .lt. sigmin) cycle
     leng = leng + 1
    endif
   enddo
  enddo
 enddo
allocate(polR(1:nr,1:leng), polH(1:nr,1:leng),polT(1:nr,1:leng),polc(1:nr,1:leng))
polR = 0.d0
polH = 0.d0
polT = 0.d0
polc = 0.d0

k = 0
 do ell = lmin, lmax
  ellp = ell + dell
  do ord = orders(ell,1),orders(ell,2)
   storderp = lbound(en(ellp)%ords,1)
   if (ell == ellp) storderp = ord
   do ordp = storderp,orders(ellp,2)

    if (existence(ell, ord) .and. existence(ellp,ordp)) then
     if (.not.(freqnu(ell)%ords(ord) .ge. freqmin .and. &
         freqnu(ellp)%ords(ordp) .ge. freqmin .and. &
         freqnu(ell)%ords(ord) .le. freqmax .and. &
         freqnu(ellp)%ords(ordp) .le. freqmax)) cycle

     freqdiff = abs(freqnu(ell)%ords(ord) - freqnu(ellp)%ords(ordp))
     if (freqdiff .gt. sigmax .or. freqdiff .lt. sigmin) cycle

     k = k + 1
     l = dble(ell)
     lp = dble(ellp)
     n = dble(ord)
     np = dble(ordp)

     ind = mapfwd(ell,ord)
     indp = mapfwd(ellp,ordp)
     conl = dble(omeg(l,0.0d0)**2)
     conlp = dble(omeg(lp,0.0d0)**2)

     if (.not. asymptotics) then

      overallcon =  8 * pi *  gamma_l(l) * gamma_l(lp) / ((4*pi)**0.5) 
  
      polR(:,k) = cmplx(0.5*(eigU(:,indp) * deigU(:,ind) - deigU(:,indp) * eigU(:,ind)) * bcoef(0.0d0,lp,dble(s),l,+1.0d0,2) &
             + 0.5*(eigV(:,indp) * deigV(:,ind) - deigV(:,indp) * eigV(:,ind)) * bcoef(1.0d0,lp,dble(s),l,+1.0d0,2))


      polR(:,k) = polR(:,k) * cmplx(rho * r**2 * overallcon) * (0.d0, 1.d0)

      polH(:,k) = cmplx((conl - conlp) * (eigU(:,ind) * eigU(:,indp) * bcoef(0.0d0,lp,dble(s),l,+1.0d0,2) &
            + eigV(:,ind) * eigV(:,indp)* bcoef(1.0d0,lp,dble(s),l,+1.0d0,2))  &
            + eigV(:,indp) * eigU(:,ind) * bcoef(1.0d0,lp,l,dble(s),+1.0d0,3)  &
            - eigV(:,ind) * eigU(:,indp) * bcoef(1.0d0,l,lp,dble(s),+1.0d0,3))

      polH(:,k) = polH(:,k) * cmplx(rho * r  * overallcon) * (0.d0, 1.d0)


      call dbyd1(tempa, aimag(polH(:,k)), nr, 1, 1)
      tempa = (drhodr*aimag(polH(:,k)) - tempa*r*stretch)/dble(s*(s+1))
      call dfftw_execute(planfwd)
      coeffs = coeffs * filter
      call dfftw_execute(planinv)
      polH(:,k) = polR(:,k)
      polR(:,k) = polR(:,k) + tempa*(0.d0,1.d0) 
 
      polT(:,k) = cmplx((eigU(:,ind) * eigV(:,indp) + eigV(:,ind) * eigU(:,indp) - eigU(:,ind) * eigU(:,indp) &
            - eigV(:,ind) * eigV(:,indp) *(conlp + conl - omeg(dble(s),0.0d0)**2)) * bcoef(1.0d0,lp,dble(s),l,-1.0d0,2))

      polT(:,k) = polT(:,k) * cmplx(rho * r)  * overallcon 

     else
      polR(:,k) = 0.d0
      !polR(:,k) = (1-2.*modulo(ell,2))*(0.d0,1.d0)*(eigU(:,ind)**2+l*(l+1)*eigV(:,ind)**2)*rho*r*l*(2*l/pi)**0.5
      polR(:,k) = (1-2.*modulo(ell,2))*(0.d0,1.d0)*(eigU(:,ind)*eigU(:,indp)-eigU(:,ind)*eigV(:,indp)-eigU(:,indp)*eigV(:,ind)+&
                                                   l*(l+1)*eigV(:,ind)*eigV(:,indp))*rho*r*l*(2*l/pi)**0.5
      polc(:,k) = -(1-2.*modulo(ell,2))*rho*c2*sqrt(2*l/pi)*& 
         (r*deigU(:,ind) - l*(l+1)*eigV(:,ind) + 2* eigU(:,ind))*(r*deigU(:,indp) - l*(l+1)*eigV(:,indp) + 2* eigU(:,indp)) 
     endif
    endif
   enddo
  enddo
 enddo

 write(llow,'(I3.3)') lmin
 write(lhigh,'(I3.3)') lmax
 write(dellst,'(I1.1)') dell
 print *,'Writing out radial points into file radius.fits'
 call writefits('radius.fits',r,nr,1,1)
 print *,'Writing out flow kernels into file '//'kernels_'//llow//'_to_'//lhigh//'_dell_'//dellst//'.fits'
 call writefits(instrument//'_kernels_'//llow//'_to_'//lhigh//'_dell_'//dellst//'.fits',aimag(polR),nr,leng,1)
 print *,'Writing out sound speed kernels into file '//'soundspeed_kernels_'//llow//'_dell_'//dellst//'_to_'//lhigh//'.fits'
 call writefits(instrument//'_soundspeed_kernels_'//llow//'_to_'//lhigh//'_dell_'//dellst//'.fits',real(polc),nr,leng,1)

print *,'Writing out metadata into file '//instrument//'_indices_'//llow//'_to_'//lhigh//'_dell_'//dellst
open(333,file=instrument//'_indices_'//llow//'_to_'//lhigh//'_dell_'//dellst,action='write',status='replace')
k = 0
do ell = lmin, lmax
 ellp = ell + dell
 do ord = orders(ell,1),orders(ell,2)
  storderp = lbound(en(ellp)%ords,1)
  if (ell == ellp) storderp = ord
  do ordp = storderp,orders(ellp,2)

   if (existence(ell,ord) .and. existence(ellp,ordp)) then
    freqdiff = abs(freqnu(ell)%ords(ord) - freqnu(ellp)%ords(ordp))
    if (freqdiff .gt. sigmax .or. freqdiff .lt. sigmin) cycle
    k=k+1
    print *,ell, ord,ellp, ordp, maxval(abs((polc(:,k)))),maxval(abs((polR(:,k))))
    if (asymptotics) write(333,*), ell,ord,ellp,ordp
    if (.not. asymptotics) write(333,*), ell,ord,ellp,ordp,s
   endif
  enddo
 enddo
enddo
close(333)
call exit()

end subroutine compute_kernels_all


!------------------


subroutine compute_kernels(l, n, lp, np, polR, polH, polT)

implicit none
integer s, ind, indp!, k
real*8 polR(nr,smin:smax), polH(nr,smin:smax), polT(nr,smin:smax), con2
real*8 conl, conlp, overallcon, l, m, n, lp, mp, np, conomeg

 ind = mapfwd(int(l),int(n))
 indp = mapfwd(int(lp),int(np))
 conl = dble(omeg(l,0.0d0)**2)
 conlp = dble(omeg(lp,0.0d0)**2)

 polR = 0.0
 polH = 0.0
 polT = 0.0

 do s=smin, smax
  overallcon = gamma_l(dble(s))  * 8 * pi  * gamma_l(l) &
             * gamma_l(lp) * (1 - 2*modulo(int(mp),2))

  !k = s - smin + 1 

  polR(:,s) = 0.5*(eigU(:,indp) * deigU(:,ind) - deigU(:,indp) * eigU(:,ind)) * bcoef(0.0d0,lp,dble(s),l,+1.0d0,2) &
            + 0.5*(eigV(:,indp) * deigV(:,ind) - deigV(:,indp) * eigV(:,ind)) * bcoef(1.0d0,lp,dble(s),l,+1.0d0,2)


  polR(:,s) = polR(:,s) * rho * r**2  * overallcon

  polH(:,s) = (conl - conlp) * (eigU(:,ind) * eigU(:,indp) * bcoef(0.0d0,lp,dble(s),l,+1.0d0,2) &
            + eigV(:,ind) * eigV(:,indp)* bcoef(1.0d0,lp,dble(s),l,+1.0d0,2))  &
            + eigV(:,indp) * eigU(:,ind) * bcoef(1.0d0,lp,l,dble(s),+1.0d0,3)  &
            - eigV(:,ind) * eigU(:,indp) * bcoef(1.0d0,l,lp,dble(s),+1.0d0,3) 

  polH(:,s) = polH(:,s) * rho * r * overallcon

  polT(:,s) = (eigU(:,ind) * eigV(:,indp) + eigV(:,ind) * eigU(:,indp) - eigU(:,ind) * eigU(:,indp) &
            - eigV(:,ind) * eigV(:,indp) *(conlp + conl - omeg(dble(s),0.0d0)**2)) * bcoef(1.0d0,lp,dble(s),l,-1.0d0,2) 

  polT(:,s) = polT(:,s) * rho * r * overallcon

 enddo
! print *,ind
! s= 25
! polH(:,1) = sqrt(rho)*eigU(:,ind) * eigV(:,indp) 
! call writefits('UV.fits',polH(:,1),nr,1,1)
! polH(:,1) = (rho)*eigU(:,ind) * eigU(:,indp) 
! call writefits('UU.fits',polH(:,1),nr,1,1)
! polH(:,2) = rho*eigV(:,ind) * eigV(:,indp)*(conlp + conl - omeg(dble(s),0.0d0)**2) 
! call writefits('VV.fits',polH(:,2),nr,1,1)
! print *,'dine'
! stop


end subroutine compute_kernels

!------------------

subroutine compute_lorentz_kernels(bpp,b00,b0p,bmp,kc,s, leng)

 implicit none
! integer*8 planfwd, planinv
 integer s,m,t, ordp, ord, ellp, ell, k, ind, indp, sigllps,leng
 !parameter (leng = (lmax-lmin+1)**2*(nord-ordmin)**2)
 complex*16, dimension(1:nr,1:leng) :: b00,bpp,bmp, kc, b0p
 real*8 l, lp, omeg0, omeg0p, omeg2, omeg2p
 real*8 b0slpl, b1slpl, b2slpl, omeg0s, omeg2s, confun, di0(1:nr)

! call dfftw_plan_dft_r2c_1d(planfwd,nr,drb0p(1), coeffs(1), FFTW_ESTIMATE)
! call dfftw_plan_dft_c2r_1d(planinv,nr,coeffs(1), drb0p(1), FFTW_ESTIMATE)

! drb0p = 0.d0
! coeffs = 0.d0
! do k=1,floor(nr/2.d0)+1
!  filter(k) = cmplx(1.d0/(1.d0 + exp(k - 1 - nr/10.d0)), 0.d0)/dble(nr)
! enddo
 
 k = 0
 do ordp = ordmin, ordmax!nord-1
  do ord = ordmin, ordmax!nord-1
   do ellp = lmin,lmax
    do ell = lmin, lmax

     k = k + 1
     l = dble(ell)
     lp = dble(ellp)

     b0slpl = bcoef(0.d0,dble(s),lp,l,1.d0,1)
     b1slpl = bcoef(1.d0,dble(s),lp,l,1.d0,1)
     b2slpl = bcoef(2.d0,dble(s),lp,l,1.d0,1)

     ind = mapfwd(ell,ord)
     indp = mapfwd(ellp,ordp)
     omeg0 = omeg(l,0.0d0)
     omeg0p = omeg(lp,0.0d0)
     omeg0s = omeg(dble(s),0.d0)

     omeg2 = omeg(l,2.0d0)
     omeg2p = omeg(lp,2.0d0)
     omeg2s = omeg(dble(s),2.d0)


     sigllps = 1-2*modulo(int(l+lp+s),2)

     confun = -4*pi*gamma_l(l)*gamma_l(lp)
              
!*(1-2*modulo(mp,2))*wigner(s,ell,ellp)%sub3j(t,m)

     bpp(:,k) = rho * confun*r*eigV(:,ind)*deigU(:,indp)*b2slpl/(4.*omeg0s*omeg2s)

     bmp(:,k) = rho * confun*r**2*deigU(:,indp)*deigU(:,ind)*b0slpl

     b0p(:,k) = rho * sigllps * confun*b1slpl/(2*omeg0s)*((r*deigU(:,ind) + &
     2.*(omeg0p**2 * eigV(:,indp) - eigU(:,indp)))*(eigU(:,ind) - eigV(:,ind)) &
      - r**2*deigU(:,indp)*deigV(:,ind))


     b00(:,k) = - rho * confun*b0slpl*&
          (r*deigU(:,ind) + 2*(omeg0**2*eigV(:,ind) - eigU(:,ind))) * &
          (r*deigU(:,indp) + 2*(omeg0p**2*eigV(:,indp) - eigU(:,indp)))

     kc(:,k) = b00(:,k)*2*gamma_l(dble(s))*c2
     
!     call dbyd1(drb0p, real(b0p/r**3), nr, 1, 1)

!     drb0p = rho * drb0p * stretch * 0.5 * r**4 /omeg0s
!     di0 = drb0p

!     call filter_vars(drb0p, di0)
!     call dfftw_execute(planfwd)
!     coeffs = coeffs*filter
!     call dfftw_execute(planinv)
!     drb0p = di0 
!     b00(:,k) = rho * b00(:,k) +  drb0p !- (3*b0p - r*drb0p)*0./omeg0s)
!     if (k==58) call writefits('divfix.fits',drb0p,nr,1,1)
!     if (k==58) call writefits('divfix0.fits',di0,nr,1,1)

    enddo
   enddo
  enddo
 enddo

! call dfftw_destroy_plan(planfwd)
! call dfftw_destroy_plan(planinv)

!stop
end subroutine compute_lorentz_kernels

!------------------

       function bcoef(N,lp,lpp,l,sig,pos)

       implicit none
       real*8 N, lp, lpp, l, sig, sig2, sig3, cg, bcoef
       integer i, pos
!       real*8, external :: factorial, wig3j
     
       bcoef = 0.d0
       sig2 = -1.0d0
       sig3 = -1.0d0
       if (modulo(int(l + lp +lpp),2)==0) sig2 = -sig2
       if (modulo(int(N),2)==0) sig3 = -sig3
       

!       cg = wig3j(lp,lpp,l,-N,0.0d0,N)

       if (pos==1) then
            cg = wigner(int(lp),int(lpp),int(l))%sub3j(-int(N), 0)
       elseif (pos==2) then
            cg = wigner(int(lpp),int(l),int(lp))%sub3j(0, int(N))
       elseif (pos==3) then
            cg = wigner(int(l),int(lp),int(lpp))%sub3j(int(N), -int(N))
       endif

       bcoef = 0.5*(1 + sig*sig2) *  sig3 * cg
       if (bcoef .ne. 0.d0) then
         do i=floor(-N +1), floor(N)
          bcoef = bcoef * ((lp+i) * (l+i))**0.5
         enddo
!         do i=floor(l -N +1), floor(l + N)
!          bcoef = bcoef * i**0.5
!         enddo
        endif


!       	( factorial(lp + N) * factorial(l+N) / &
!     	&   (factorial(lp - N) * factorial(l-N)))**0.5
        

      end function bcoef
      
!------------------


      function factorial(N)

      implicit none
      real*8 factorial, N
      integer i

      factorial = 1.0d0
      if (N > 1) then
       do i=2,int(N)
        factorial = factorial * i
       enddo
      endif
      
      end function factorial 

!------------------

      function gamma_l(ell)

      implicit none
      real*8 gamma_l, N, ell
      integer i

      gamma_l = ((2*ell+1)/(4*pi))**0.5
      
      end function gamma_l 


!------------------


      function omeg(l, N)
      implicit none
      real*8 omeg, l, N
       omeg = (0.5*(l+N)*(l-N+1))**0.5
      end function omeg

!------------------

subroutine BINARYREADER

 implicit none
! include 'params.i'
 real*4 Grav, arr(100), l, zk, muhz
 real*4 ell, freq
 real*4, dimension(:,:), allocatable :: X
 real*8 arr2(100), denav, norm
 real*8, dimension(:,:), allocatable :: vars
 integer i, n, norder

 ordmin = 1
 ordmax = 10

 allocate(X(4,nr), vars(6,nr))

 open(99, file='/scratch/jb6888/QDP/sfopal5h5', form='unformatted', status='old')
 read(99)
 read(99) arr2

 rsun= arr2(1)
 msun = arr2(2)
 rhonorm = msun/(4.*dble(pi)*rsun*rsun*rsun)

 freqnorm = sqrt(Gcon*msun/(rsun*rsun*rsun))
 speednorm = sqrt(Gcon*msun/rsun)

 do i=1,nr
  read(99) vars(:,i)
   
  r(i) = vars(1,i)
  rho(i) = vars(2,i)*vars(6,i)
  c2(i) = r(i)**2*vars(2,i)/vars(3,i)
 enddo
 close(99)
 r = r(nr:1:-1)
 dr(1:nr-1) = r(2:nr) - r(1:nr-1)
 dr(nr) = dr(nr-1)

 rho = rho(nr:1:-1)
 c2 = c2(nr:1:-1)

 open(99, file='/scratch/jb6888/QDP/egvt.sfopal5h5', form='unformatted', status='old')
 read(99)
 read(99) arr

 n = 0
 l = 0
 mapfwd = 0
 do while (int(l) .le. lmax)
  read(99)  norder, l, zk, muhz
  read(99) X
  if ((int(l) .ge. lmin) .and.  (int(l) .le. lmax)) then
   if ((norder .le. orders(int(l),2)) .and. (norder .ge. orders(int(l),1))) then
    if (existence(int(l),norder)) then
      n = n + 1
      eigU(:,n) = dble(X(1,nr:1:-1))
      eigV(:,n) = (dble(X(2,nr:1:-1))/rho + dble(X(3,nr:1:-1)))/(r*zk**2) !*(1.-r)**2!/dble(l*(l+1))
      norm = sum(dr*r**2*rho*(eigU(:,n)**2 + l*(l+1)*eigV(:,n)**2))**0.5
      eigU(:,n) = eigU(:,n)/norm
      eigV(:,n) = eigV(:,n)/norm
      mapfwd(int(l), norder) = n
       !if (int(l)==30) print *,norder,l,zk,muhz
      !omegaref(n) = muhz * 1e-6 * 2.0 * pi!/freqnorm
      !print *,omegaref(n),muhz,minval(r),maxval(r)
    endif
   endif
  endif
 enddo
 
! do i=1,nr
!  print *,r(i),c2(i)**0.5*speednorm*1e-5/695989.4
! enddo
! stop
 deallocate(X, vars)
END subroutine BINARYREADER

! ==========================================

        SUBROUTINE writefits(filename,temp,nx, ny, nz)

        implicit none 
        integer blocksize,bitpix,naxes(3),unit1,nz , nx, ny
        integer status1,group,fpixel,flag, nelements
	  real*8 temp(nx,ny,nz)
!	  real*8 dump_array(nx,ny,nz)
	  character*(*) filename
	  logical simple,extend, exists

          inquire(file=filename, exist=exists)
          if (exists) call system('rm '//filename)
          print *,'Writing into file '//filename

	  status1 = 0
	  call ftgiou(unit1,status1)
	  blocksize=1
!	 dump_array = dble(temp)
	  call ftinit(unit1,filename,blocksize,status1)
	  simple=.true.
	  bitpix=-64
	  naxes(1)=nx
	  naxes(2)=ny
	  naxes(3)=nz
	  nelements=naxes(1)*naxes(2)*naxes(3)
	  extend=.false.
	  group=1
	  fpixel=1

	  call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
	  call ftpprd(unit1,group,fpixel,nelements,temp,status1)
	  call ftclos(unit1, status1)
	  call ftfiou(unit1, status1)


	 end SUBROUTINE writefits
!================================================================================

	 SUBROUTINE readfits(filename,readarr,dim1,dime2,dim3)

	  implicit none
	  integer status,unit,readwrite,blocksize,naxes(3)
	  integer group,firstpix, dim3 ,dime2,dim1
	  integer nelements
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
	 function sting(numbe,length_string)

	  implicit none
	  integer i,length_string,numbe,n(1:length_string),number_temp
	  character*(length_string) sting
	  character*1 charc(10)

	  charc(1)  = '0'
	  charc(2)  = '1'
	  charc(3)  = '2'
	  charc(4)  = '3'
	  charc(5)  = '4'
	  charc(6)  = '5'
	  charc(7)  = '6'
	  charc(8)  = '7'
	  charc(9)  = '8'
	  charc(10) = '9'


	  number_temp = numbe
	  do i=length_string,1,-1

	   n(length_string-i+1) = floor(number_temp*10.0**(-(i-1.0)))
	   number_temp = number_temp - n(length_string-i+1)*10**(i-1)
 	   sting(length_string-i+1:length_string-i+1)  &
              = charc(n(length_string-i+1)+1)

	  enddo

	  end function sting

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
SUBROUTINE FILTER_VARS(inp,outp)
  implicit none
  real*8 coeffs(0:5), inp(nr), outp(nr)
  integer st, fi, l, rad
                                                                                                                                
  coeffs(0) = 0.75390625000000
  coeffs(1) = 0.41015625000000
  coeffs(2) = -0.23437500000000
  coeffs(3) = 0.08789062500000
  coeffs(4) = -0.01953125000000
  coeffs(5) = 0.00195312500000
         
  outp = 0.d0
  Rad = nr                                                                                                                       
  outp(6:Rad-5) = coeffs(0)*inp(6:Rad-5)  &
     +  0.5*coeffs(1)*(inp(5:Rad-6) + inp(7:Rad-4)) &
     +  0.5*coeffs(2)*(inp(4:Rad-7) + inp(8:Rad-3))  &
     +  0.5*coeffs(3)*(inp(3:Rad-8) + inp(9:Rad-2))  &
     +  0.5*coeffs(4)*(inp(2:Rad-9) + inp(10:Rad-1))  &
     +  0.5*coeffs(5)*(inp(1:Rad-10) +inp(11:Rad))

END SUBROUTINE FILTER_VARS

!================================================================================================
subroutine read_leakage(dl,dm)
  implicit none
  integer dm, dl, dl_mat,dm_mat, ind, ind0, i, j, l, lp, m, mp
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

   call readfits('/scratch/jb6888/QDP/leakvw0/default/vradsum/leakrlist.vradsum.fits',cr_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)
   call readfits('/scratch/jb6888/QDP/leakvw0/default/vradsum/leakilist.vradsum.fits',ci_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)

   call readfits('/scratch/jb6888/QDP/leakvw0/default/vhorsum/leakrlist.vhorsum.fits',hr_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)
   call readfits('/scratch/jb6888/QDP/leakvw0/default/vhorsum/leakilist.vhorsum.fits',hi_mat,&
      2*dm_mat+1,2*dl_mat+1,ind)

   
   allocate(rleaks(-dm:dm,-dl:dl,lmin:lmax), horleaks(-dm:dm,-dl:dl,lmin:lmax))
   do lp = lmin, lmax
    ind0 = lp*(lp+1)/2
    do j = -dl, dl
     do i = -dm, dm
      allocate(rleaks(i,j,lp)%em(-lp:lp),horleaks(i,j,lp)%em(-lp:lp))
     enddo
    enddo

    do mp = 0, lp
     ind = ind0 + mp + 1
     do j = -dl, dl
      l = lp + j

      do i = -dm, dm

       m = mp + i
       if (abs(m) > l) cycle

       L_ii = 0.5*sign(1,m)*ci_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)
       L_rr = 0.5*cr_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)

       rleaks(i,j,lp)%em(mp) = L_rr + L_ii  
       rleaks(i,j,lp)%em(-mp) = L_rr + sign(1,mp)* L_ii

       L_ii = 0.5*sign(1,m)*hi_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)
       L_rr = 0.5*hr_mat(abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind)

       horleaks(i,j,lp)%em(mp) = L_rr + L_ii  
       horleaks(i,j,lp)%em(-mp) = L_rr + sign(1,mp)* L_ii

      enddo
     enddo
    enddo
   enddo

   deallocate(cr_mat,ci_mat,hr_mat,hi_mat)

   !leak = (cr+ci)/2.d0

!  else

   !leak = 0.d0
!   print *,'BEYOND ACCEPTABLE RANGE!'

!  endif

 

end subroutine read_leakage

!================================================================================================

function dfactorial(q)
  implicit none
  integer q, i
  real*8 dfactorial

  dfactorial = 1.d0
  do i=q,1,-2
   dfactorial = dfactorial * i
  enddo

end function dfactorial

!================================================================================================

subroutine compute_Q3(dm, dl, dtee, des, tee, ess,Qfun,sig, year,cond)

 implicit none

 integer dl, es, ellp, ell, j, ems, emsp, emsdp, rind, sp, dell, dtee
 integer elldp, dm, delt, i, h, esp, t, tp, slow, shigh, tee, ess, ierr, cond
 integer tmin,tmax, IER, des, k, l, s, q, ord, ind, lbo, ubo, year
 integer intarr(1:1500), delta_nu, sigoff, splow, sphigh
 real*8 lmi, lma, temv(1:(lmax+lmax+1)), func(1:nr), f_sq, leak1, mix1
 real*8 leak2, mix2, constantwopi, twopiemin6
 real*8 fofr(1:nr), omegal, con1, alphanla, nyears
 real*8, dimension(:), allocatable, target :: nu, mask
 real*8, dimension(:,:), allocatable :: norms
 real*8, dimension(:), pointer :: nuf
 complex*16, dimension(:), allocatable, target :: array1, array2, array3, powpos, powneg
 complex*16, dimension(:), pointer :: powerfreq, powerfreq_p, hnl
 real*8 nj(-dl:dl), gamnl, gamnlp, nulmn,nulmnt, contwopi, dnu, conabs(-dtee:dtee,-des:des)
 real*8 lolim_mp,hilim_mp, lolim_m, hilim_m, Csum, ref, norm, sig, con0
 integer  minoff, maxoff!, dl, dm
 complex*16 con(-dtee:dtee,-des:des), compar, hnlp,Qfun(lmin+6:lmax-6,0:30,1:(2*des+1)*(2*dtee+1))
 character*2 tc, sc, spc, tpc,ynum
 character*3 lch_p
 logical lexist

 t = tee
 s = ess
 sphigh = min(s+des,smax)
 contwopi = 1.d0/(2.d0*pi)**2
 constantwopi = 1e12 * contwopi
 twopiemin6 = 2*pi*1e-6

 write(ynum,'(I2.2)') year

 nyears = 1.0
 allocate(gammast(lmin:lmax), gammasptp(-2*dl:2*dl,lmin:lmax,-des:des))
  
 do ellp = lmin,lmax

  allocate(gammast(ellp)%sub3j(-ellp:ellp))
  gammast(ellp)%sub3j(:) = 0.d0

  do ems = -ellp, ellp

   temv = 0.d0
   t = tee
   if (abs(ems+t) .gt. ellp) cycle
   call drc3jj(dble(ellp),dble(ellp),dble(ems),dble(-ems-t),lmi,lma,temv,lmax+lmax+1,IER) 

   if (ess - int(lmi) + 1 > 0) &
     gammast(ellp)%sub3j(ems) = temv(ess-int(lmi) + 1) * &
                                      (2.d0*ess+1.d0)**0.5*(1.-2.*modulo(tee+ems,2)) 

  enddo

  do ell = max(lmin,ellp-2*dl),min(lmax,ellp+ 2*dl)
   j = ell - ellp

   do sp = max(s - des,smin), sphigh
    tmin = max(tee-dtee,-sp)
    tmax = min(tee+dtee,sp)

    allocate(gammasptp(j,ellp,sp-s)%sub3j(tmin:tmax,-ell:ell))
    gammasptp(j,ellp,sp-s)%sub3j(tmin:tmax,:) = 0.d0
    do tp = tmin, tmax

     do ems = -ell, ell

      if (abs(ems+tp) .gt. ellp) cycle
      temv = 0.d0
      call drc3jj(dble(ell),dble(ellp),dble(ems),dble(-ems-tp),lmi,lma,temv,lmax+lmax+1,IER) 

      if (sp - int(lmi) + 1 > 0) &
       gammasptp(j,ellp,sp-s)%sub3j(tp,ems) = temv(sp-int(lmi) + 1) * &
                                      (2.d0*sp+1.d0)**0.5*(1.-2.*modulo(tp+ems,2))
     enddo
    enddo
   enddo

  enddo
 enddo
 s = ess
 t = tee
 esp = sp

 Qfun = 0.d0


 dnu = 1./(86400.*360.) * 1e6
 norm = 2.5e-10 * 2.0 / nyears

 sigoff = floor((sig + t * (trackrate-0.0317))/dnu)!a1(ell)%ords(ord))/dnu) !

 minoff = floor((freqmin-150)/dnu)
 maxoff = floor((freqmax+150)/dnu)
 allocate(nu(minoff:maxoff),array1(minoff:maxoff),mask(minoff:maxoff),&
         array2(minoff:maxoff), array3(minoff:maxoff),powpos(minoff:maxoff),&
         powneg(minoff:maxoff))

 nuf => nu(minoff:maxoff)
 powerfreq => array1(minoff:maxoff)
 powerfreq_p => array2(minoff:maxoff)
 hnl => array3(minoff:maxoff)

 if (t > 0) then
  allocate(leakmat(-dm:dm,-dl:dl))
  do dell = -dl, dl
   do i = -dm,dm
    allocate(leakmat(i,dell)%em(-(lmax-6):(lmax-6)))
   enddo
  enddo
 endif


 Qfun = 0.d0
 nu = (/(i*dnu, i=minoff,maxoff)/)

 conabs = 0.d0

  
 do l = lmin+6, lmax-6
  print *,l

 
  if (allocated(norms)) deallocate(norms)
  allocate(norms(lbound(en(l)%ords,1):ubound(en(l)%ords,1),-dl:dl))
  norms = 0.d0!1.35e-10/nyears

  do dell = -dl,dl
   elldp = l + dell
   
   write(lch_p,'(I3.3)') elldp
   inquire(file='/scratch/shravan/'//instrument//'/norms/'//instrument//&
           '_'//lch_p//'_year_'//ynum,exist=lexist)
   if (lexist) then
    open(1555,file='/scratch/shravan/'//instrument//'/norms/'//instrument//&
          '_'//lch_p//'_year_'//ynum,status='old',action='read')
    do  
     read(1555,*,IOSTAT=ierr) con0, ord
     if (ierr .ne.0) exit
   
     if (ord .ge. lbound(en(l)%ords,1) .and. ord .le. ubound(en(l)%ords,1) &
     .and. con0 .lt. 1e-6 .and. (.not. isnan(con0))) &
      norms(ord,dell) = con0/nyears
    enddo
    close(1555)
   endif
  enddo
!cycle
  do ord = lbound(en(l)%ords,1),ubound(en(l)%ords,1)
   if (.not.existence(l,ord)) cycle
   do dell = -dl, dl
    mix1 = 0.d0
    elldp = l + dell
    nj(dell) = 0.d0
    do i=-dm,dm
     leakmat(i,dell)%em = 0.d0
    enddo
    if (existence(elldp,ord)) then
     mix1 = (274.8*1e2*(elldp+0.5)/((twopiemin6)**2*rsun))/freqnu(elldp)%ords(ord)**2
     nj(dell) =norms(ord,elldp-l) * fwhm(elldp)%ords(ord) * (amps(elldp)%ords(ord) * freqnu(elldp)%ords(ord))**2 * (2*pi)**3 * 1e-18
    endif

    do i = -dm,dm
     leakmat(i,dell)%em(-l:l) = rleaks(i,dell,l)%em(-l:l) + mix1 * horleaks(i,dell,l)%em(-l:l)
    enddo

   enddo

   Csum = 0.d0
   con = (0.d0, 0.d0)
   gamnl = fwhm(l)%ords(ord)

   !if (nus(l)%mn(-l,ord)-gamnl-sig < freqmin .or. nus(l)%mn(l,ord)+gamnl > freqmax  .or. 
   if (norms(ord,0)==0.d0) cycle
   do ems=max(-l-t,-l),min(l,l-t)
    nulmn = nus(l)%mn(ems,ord)
    nulmnt = nus(l)%mn(ems+t,ord)
    lolim_m = nulmn - gamnl
    hilim_m = nulmn + gamnl
    lolim_mp = nulmnt - gamnl - (sig + t * (trackrate-0.0317))
    hilim_mp = nulmnt + gamnl - (sig + t * (trackrate-0.0317))
      
    lolim_m = min(lolim_m,lolim_mp)
    hilim_m = max(hilim_m,hilim_mp)
    minoff = floor(lolim_m/dnu)-1
    maxoff = floor(hilim_m/dnu)+1 

    powpos = 1.d0/(nu + nulmn)
    powneg = 1.d0/(nu - nulmn + (0.d0,0.5d0)*gamnl)

    delta_nu = floor((-nulmn+nulmnt)/dnu)

    mask(minoff:maxoff) = 0.d0
    where ((((lolim_mp .le. nu(minoff:maxoff)) .and. (hilim_mp .ge. nu(minoff:maxoff))) &
    .or. ((lolim_m .le. nu(minoff:maxoff)) .and. (hilim_m .ge. nu(minoff:maxoff))))) &
                mask(minoff:maxoff) = 1.d0

    k = 0
    do j=minoff,maxoff
     if (mask(j) == 1.d0) then
      k = k+1
      intarr(k) = j
     endif
    enddo
    if (k==0) cycle

     
    powerfreq(1:k) = - powpos(intarr(1:k)) * powneg(intarr(1:k)) * 1e12 * contwopi

    powerfreq_p(1:k) = -powpos(intarr(1:k) + delta_nu + sigoff) * powneg(intarr(1:k) - delta_nu + sigoff) * 1e12 * contwopi 

    leak1 = leakmat(0,0)%em(ems)  * leakmat(0,0)%em(ems+t)  
    hnl(1:k) = -2*nu(intarr(1:k))*twopiemin6*nj(0)*conjg(powerfreq_p(1:k)*abs(powerfreq(1:k))**2 + &!
                     conjg(powerfreq(1:k))*abs(powerfreq_p(1:k))**2) * leak1 * gammast(l)%sub3j(ems)

    Csum = sum(abs(hnl(1:k))**2) + Csum

    do  elldp = max(lmin,l-dl), min(lmax,l+dl)

     if (lbound(nus(elldp)%mn,2) > ord .or. ubound(nus(elldp)%mn,2) < ord) cycle
     gamnlp = fwhm(elldp)%ords(ord)
     
     !if (nus(elldp)%mn(-elldp,ord)-gamnlp-sig < freqmin .or. nus(elldp)%mn(elldp,ord)+gamnlp > freqmax .or. &
      if ( norms(ord,elldp-l)==0.d0 ) cycle

      do ellp = max(lmin,l-dl), min(lmax,l+dl)
       if (lbound(nus(ellp)%mn,2) > ord .or. ubound(nus(ellp)%mn,2) < ord) cycle 

       gamnl = fwhm(ellp)%ords(ord)
       q = elldp - ellp
       !if (nus(ellp)%mn(-ellp,ord)-gamnl-sig < freqmin .or. nus(ellp)%mn(ellp,ord)+gamnl > freqmax &
        if( (abs(q) > dl)  .or. norms(ord,ellp-l)==0.d0  ) cycle


       do tp = t-dtee, t+dtee
        if (tp==0) cycle

        lbo = max(ems-dm,-ellp,-tp-elldp,-dm+ems+t-tp)
        ubo = min(ems+dm,ellp,elldp-tp,dm+ems+t-tp)

        do emsp=lbo,ubo
         emsdp = (emsp + tp)

         delta_nu = floor((nus(ellp)%mn(emsp,ord) - nulmn)/dnu)
         powerfreq(1:k) = - powpos(intarr(1:k) + delta_nu) * powneg(intarr(1:k) - delta_nu) * constantwopi

         delta_nu = floor((nus(elldp)%mn(emsdp,ord) - nulmn)/dnu)
         powerfreq_p(1:k) = -powpos(intarr(1:k) + delta_nu + sigoff) * powneg(intarr(1:k) -delta_nu +sigoff) * constantwopi

         leak2 = leakmat(emsdp-ems-t,elldp-l)%em(ems+t) *  leakmat(emsp-ems,ellp-l)%em(ems) 

         hnlp = sum(-2*nu(intarr(1:k))*twopiemin6 *(nj(ellp-l)*powerfreq_p(1:k)*abs(powerfreq(1:k))**2 + & !
                     nj(elldp-l)*conjg(powerfreq(1:k))*abs(powerfreq_p(1:k))**2)  *hnl(1:k)) * &
                      leak2
         splow = max(s-des,smin,abs(tp))
         do sp = splow, sphigh
          if (modulo(sp + ellp + elldp,2).ne.cond) cycle
          f_sq = dfactorial(sp - q+cond-1) * dfactorial(sp + q+cond-1) / (sqrt(factorial(dble(sp-q))*factorial(dble(sp+q)))) &
           *(1 - 2*modulo((sp+elldp-ellp-cond)/2,2))
          con(tp-t,sp-s) = gammasptp(ellp-elldp,elldp,sp-s)%sub3j(tp,emsp) *  hnlp * f_sq + con(tp-t,sp-s)
         enddo !sp
        enddo !tp
       enddo ! emsp

      enddo ! ellp
     enddo ! elldp
    enddo ! ems

    do sp = s-des,s+des
     do tp = t-dtee,t+dtee
      Qfun(l,ord,(2*dtee+1)*(sp-s+des) + tp - t + dtee + 1) =  con(tp-t,sp-s)/Csum
     enddo
    enddo
    conabs = conabs + abs(con/Csum)**2

   enddo ! ord
  enddo ! l

 if (1==2) then
  write(sc,'(I2.2)') ess 
  write(tc,'(I2.2)') abs(tee) 
  write(spc,'(I2.2)') sp
  write(tpc,'(I2.2)') abs(tp)
  if (tee > 0)  open(233,file='/scratch/shravan/'//instrument//'/kernels/relpows_'//&
           tc//'_'//sc,status='unknown',position='append') 

  if (tee < 0)  open(233,file='/scratch/shravan/'//instrument//'/kernels/relpows_neg_'//&
           tc//'_'//sc,status='unknown',position='append') 


  do sp = s-des,s+des
   do tp = max(t - dtee,-sp),min(t+dtee,sp)
    write(233, *)conabs(tp-t,sp-s),ess,tee,sp,tp,(ess-sp),(tee-tp)
    print *,conabs(tp-t,sp-s),ess,tee,sp,tp,(ess-sp),(tee-tp)
   enddo
  enddo
  flush(233)
  close(233)

 endif

  deallocate(gammast, gammasptp,array1, array2, array3, nu)
  if (t < 0) then
   do dell = -dl, dl
    do i = -dm,dm
     deallocate(leakmat(i,dell)%em)
    enddo
   enddo
  endif

 end subroutine compute_Q3

!================================================================================================
 end module kernels
