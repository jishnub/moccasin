 PROGRAM INVERSION

  use kernels
  use bspline
  use frequencies

  implicit none
  integer ell, s, ord, leng, one, nrad, info, i, ellp, ordp, nyears,sp,tp, mdp
  integer kord, lwork, k, ind, sstart, indlow,indspace, t, sigoff,m,mp, elldp
!  integer ordmin, ordmax, nord
  integer, allocatable, dimension(:) :: ipiv
  real*8 sig, offset, rlower!, pow(smin:smax, ordmin:ordmax)
  real*8 es, l, lp, tstart, tfin!, tempar(609,305)
  complex*16 zeroc, alph, reg
  real*8, allocatable, dimension(:) :: rknot, rsmall, &
                                      drsmall, temp
  complex*16, allocatable, dimension(:) :: rhs, rhsout, &
                                    work, output_spline
  complex*16,allocatable, dimension(:,:) :: polH, polT, polR, kc, k0p
  complex*16,allocatable, dimension(:,:) ::  Amat !, Aout
  real*8, allocatable, dimension(:,:) :: aout
  real*8, allocatable, dimension(:) :: qfun
  complex*16, allocatable, dimension(:,:,:) :: tempb, Lfun
  character*1 charc, charn
  character*2 str, yearc,tpc,spc,tc,sc
  character*3 ellc, trackch
  integer ilaenv, yearnum
  external ilaenv
  parameter(one = 1, zeroc = (0.d0, 0.d0), charc = 'C')
  parameter(alph = (1.d0, 0.d0), charn = 'N', kord = 5, nrad=100) 
  logical sound_speed_kernels
  

  real*8, dimension(:), allocatable :: nu
  complex*16, dimension(:), allocatable :: powerfreq, powerfreq_p, hnl, hnlp
  real*8 nj, njp, gamnl, gamnlp, nulmn,nulmnt, con, contwopi, dnu
  real*8 lolim_mp,hilim_mp, lolim_m, hilim_m
  integer  minoff, maxoff, dl, dm, cond



  ell=60
  m=10
  mp = 32
  sound_speed_kernels = .false.

  write(trackch,'(I3.3)') nint(1e3*trackrate)
  cond = 1
  if(sound_speed_kernels) cond = 0
  print *,'Standard basic setup'
  call basic_setup
  !MODE_ANALYSIS = .FALSE.
 if (1==2) then

  call getarg(1,tc)
  read(tc,*) t

  call getarg(2,sc)
  read(sc,*) s

  print *,'Finished setup'
!  call construct_frequencies('/scratch/data/4816/m10q.36',lmin,lmax,3)

  call read_leakage(6,15)
!  ell = 126
!  ellp=125
!  elldp = 125

!  dl = 6
! con = 2.5e-10
! contwopi  = 1.d0/(2*pi)**2
! ord = 2
! sig = 4.d0
!  m = 10
!  mp = 20
!  dm = 15


  !call compute_Q(15,6,10,10)
  allocate(qfun(1:nr),Lfun(lmin+6:lmax-6,0:26,1:63))

!  es=0.d0
!  s = 21
!  t= 12
  write(sc,'(I2.2)') s 
  write(tc,'(I2.2)') t 
   
!  print *,'here'
  Lfun = 0.d0
  print *,t,s
  yearnum = 1
  call compute_Q3(4,2,4,3,t,s,Lfun,1.d0,yearnum,cond)
!  call compute_Q3(4,2,2,2,t,s,Lfun,1.d0,yearnum)
  

 call writefits('/scratch/shravan/'//instrument//'/kernels/tracking'//trackch//'/Q_real_'//& !_sound
                     tc//'_'//sc//'.fits',real(Lfun),lmax-lmin-11,27,63)

 call writefits('/scratch/shravan/'//instrument//'/kernels/tracking'//trackch//'/Q_imag_'//& !_sound
                     tc//'_'//sc//'.fits',aimag(Lfun),lmax-lmin-11,27,63)

  t = -t
  call compute_Q3(4,2,4,3,t,s,Lfun,1.d0,yearnum,cond)
  

 call writefits('/scratch/shravan/'//instrument//'/kernels/tracking'//trackch//'/Q_real_neg_'//& !_sound
                     tc//'_'//sc//'.fits',real(Lfun),lmax-lmin-11,27,63)

 call writefits('/scratch/shravan/'//instrument//'/kernels/tracking'//trackch//'/Q_imag_neg_'//& !_sound
                     tc//'_'//sc//'.fits',aimag(Lfun),lmax-lmin-11,27,63)

 call system('rm /scratch/shravan/kernprocess/'//tc//'_'//sc)

  call exit()
  stop
 endif

 
if (1==1) then

!  ordmax = 15
!  ordmin = 1
  nord = ordmax - ordmin + 1
  leng = ntot!sum(orders(:,2) - orders(:,1)) + lmax - lmin + 1!(lmax-lmin+1)*(ordmax-ordmin+1)
  allocate(polR(1:nr,1:leng), kc(1:nr,1:leng), &
          polH(1:nr,1:leng), polT(1:nr,1:leng), &
          k0p(1:nr,1:leng))
  es=29.d0
!  call compute_lorentz_kernels(polT, polR, k0p, polH, kc, int(es))
 call compute_kernels_inversion(polR, polH, polT, int(es), leng)

  ! print *,maxval(abs(polR)),maxval(abs(polH)),maxval(abs(polT))
!  print *,leng,nr,maxval(abs(polR)),maxval(abs(polH))

! stop
!  print *,'here', leng
  call writefits('radius.fits',r,nr,1,1)
  call writefits(instrument//'_kernels.fits',aimag(polR),nr,leng,1)
  call writefits(instrument//'_soundspeed_kernels.fits',real(polT),nr,leng,1)
!  call writefits(instrument//'_true_kernels.fits',real(polT),nr,leng,1)
  !call writefits('kernR_elastic_ell_120_121_es21_n0to6.fits',aimag(polH),nr,leng,1)
!  call writefits('kernT_ell_120_121_es21_n0to6.fits',real(polT),nr,leng,1)
 
 if (1==2) then
  call writefits('kern00_ell_60_61_es21_n0to6.fits',real(polR),nr,leng,1)
  call writefits('kern0+_ell_60_61_es21_n0to6.fits',real(k0p),nr,leng,1)
  call writefits('kern+-_ell_60_61_es21_n0to6.fits',real(polH),nr,leng,1)
  call writefits('kern++_ell_60_61_es21_n0to6.fits',real(polT),nr,leng,1)
  call writefits('kernc_ell_60_61_es21_n0to6.fits',real(kc),nr,leng,1)

 endif
if (1==2) then
  call writefits('kern00_ell_120_121_es21_n0to6.fits',real(polR),nr,leng,1)
  call writefits('kern0+_ell_120_121_es21_n0to6.fits',real(k0p),nr,leng,1)
  call writefits('kern+-_ell_120_121_es21_n0to6.fits',real(polH),nr,leng,1)
  call writefits('kern++_ell_120_121_es21_n0to6.fits',real(polT),nr,leng,1)
  call writefits('kernc_ell_120_121_es21_n0to6.fits',real(kc),nr,leng,1)
endif

 if (1==2) then
do s = 4,35
  print *,'Doing s:',s
  !call compute_kernels_inversion(polR, polH, polT, s)
  write(str,'(i2)') s
  if (s<10) str(1:1) = '0'
  call writefits('kern_U_ell_100_120_es_'//str//'_n_2.fits',aimag(polR),nr,leng,1)
  call writefits('kern_W_ell_100_120_es_'//str//'_n_2.fits',real(polT),nr,leng,1)
  call writefits('asymp_ell_100_120_es_'//str//'_n_2.fits',aimag(polH),nr,leng,1)
enddo


endif
open(333,file=instrument//'_indices',action='write',status='replace')
s = 0
!do ordp = ordmin, ordmax!nord-1
!  do ellp = lmin,lmax
 do ell = lmin+6, lmax-6
  do ord = orders(ell,1),orders(ell,2)
   if (existence(ell,ord)) then
    ellp = ell
    s=s+1
    print *,ellp, ell, ord, ordp, s,maxval(abs((polT(:,s)))),maxval(abs((polR(:,s))))!,&
              !maxval(abs((polH(:,s)))),maxval(abs((kc(:,s)))),maxval(abs(k0p(:,s)))
    write(333,*), ellp,ell,ord,ordp,s
    endif
   enddo
!  enddo
 enddo
!enddo
close(333)

  stop

  offset = -0.5d0

  l = 120.d0
  lp = 120.d0
  !print *,'here' ,2.d0*omeg(l,0.d0)**2,wigner(int(es),int(lp))%sub3j(1,1)
  print *,bcoef(1.d0,l,es,lp,-1.d0,2),bcoef(1.d0,es,l,lp,-1.d0,1) 
  print *,
  print *, omeg(l,0.d0)**2 * bcoef(1.d0,es,l,lp,-1.d0,1) + &
 2.*modulo(int(es+l+lp),2)*omeg(l,0.d0)*omeg(lp,0.d0)*omeg(l,2.d0)*omeg(es,0.d0)*wigner(int(es),int(l),int(lp))%sub3j(1,-2),&
 (2*omeg(l,0.d0)**2 - omeg(es,0.d0)**2)*bcoef(1.d0,l,es,lp,-1.d0,2)
  stop

endif

  ell = 120
  sig = 2.d0
  ord = 2
  !call analyzethis (ell)
  allocate(tempb(0:32,0:smax,sigmin:sigmax))
  tempb = 0.d0
  do sigoff= sigmin,sigmax
  do t = 5,32
   do s = t, t+3
!    tempb(t,s,sigoff) = b_coefs(s,ord,120,sigoff)%tee(t)
   enddo
  enddo
!    t=7
!    s=10
!    print *,tempb(t,s,sigoff)
  enddo
  print *,'here'
  call writefits('realb.fits',real(tempb),33,smax+1,sigmax-sigmin+1)
  call writefits('imagb.fits',aimag(tempb),33,smax+1,sigmax-sigmin+1)
  deallocate(tempb)
  stop
  !call compute_power_s_n(pow)

!  open(498,file='totpower',action='write')
  do s=smin,smax
!    write(498,*) pow(s,:)
   !print *,sum(abs(b_coefs(-smax:smax,s,0:nord-1,lmin:lmax))**2)
!   write(498,*) sum(abs(b_coefs(s,0,lmin:lmax)%tee(-s:s))**2/(2*s+1))**0.5,&
!          &sum(abs(b_coefs(s,1,lmin:lmax)%tee(-s:s))**2/(2*s+1))**0.5,&
!              sum(abs(b_coefs(s,2,lmin:lmax)%tee(-s:s))**2/(2*s+1))**0.5,&
! &sum(abs(b_coefs(s,3,lmin:lmax)%tee(-s:s))**2/(2*s+1))**0.5,&
!  sum(abs(b_coefs(s,4,lmin:lmax)%tee(-s:s))**2/(2*s+1))**0.5,&
!sum(abs(b_coefs(s,5,lmin:lmax)%tee(-s:s))**2/(2*s+1))**0.5
!,sum(abs(b_coefs(-smax:smax,s,6,lmin:lmax))**2)!,&
!sum(abs(b_coefs(-smax:smax,s,7,lmin:lmax))**2)
  enddo
 ! close(49819)


  reg = (1e-8, 0.d0)

  allocate(rsmall(1:nrad), rknot(nrad+kord), drsmall(1:nrad))

  rlower = 0.9
  indlow = minloc(abs(r-rlower),dim=1)
  indspace = floor(dble(nr - indlow+1)/dble(nrad))  
  indlow = nr - indspace * (nrad - 1)
  rsmall = r(indlow:nr:indspace)
  drsmall(1:nrad-1) = r(2:nrad) - r(1:nrad-1)
  drsmall(nrad) = drsmall(nrad-1)

  call dbsnak(nrad,rsmall,kord,rknot)
  
  
  lwork = nrad*ilaenv(1, 'zgetri','n',nrad,nrad,leng,leng)

  leng = (lmax-lmin+1)*(nord-ordmin)
  allocate(Aout(1:nrad,1:nrad),Amat(1:leng,1:nrad),&
          rhs(1:leng), rhsout(1:nrad),w_st(1:nrad,smin:smax),&
          ipiv(1:nrad),work(1:lwork), polR(1:nr,1:leng), &
          polH(1:nr,1:leng), polT(1:nr,1:leng), temp(1:nrad), &
          output_spline(1:nrad))
  Amat = (0.d0, 0.d0)
 
  sstart=smin + 1 - modulo(smin,2)
 
  do s = sstart, smax, 2

   !call compute_kernels_inversion(polR, polH, polT,s)
!   Amat = polT ! if nrad = nr

   do i= 1,leng
    call dbsint(nrad,rsmall,real(polT(indlow:nr:indspace,i)*drsmall), &
              kord,rknot,temp)
    Amat(i,:) = cmplx(1.d0, 0.d0) * abs(temp)**2
     
   enddo

   print *,'finished kernels'

   Aout = (0.d0, 0.d0)
   do k=1,nrad
    Aout(k,k) = (1.d0,0.d0)
   enddo

   call zgemm(charc, charn, nrad, nrad, leng, alph, Amat, leng, &
              Amat, leng, reg, Aout, nrad)


!   print *,'finished matrix multiplication. now doing inversion'
   call zgetrf(nrad, nrad, Aout(1,1), nrad, ipiv(1), info)
!   print *,'Determined the pivots'
   call zgetri(nrad, Aout(1,1), nrad, ipiv(1), work(1), lwork, info)
!   print *,'Finished inversion'

!   do t = 1,s

    print *,'Doing s', s
    rhsout = zeroc
    rhs = zeroc

    ind = 0
    do ell=lmin,lmax
     do ord=ordmin,ordmax
      ind = ind+1
 ! NEEDS TO BE UNCOMMENTED
!      rhs(ind) = sum(abs(b_coefs(s,ord,ell)%tee(-s:s))**2)/(2.d0*s + 1)
     enddo
    enddo

    call zgemv(charc,leng,nrad,alph,Amat,leng,rhs,one,zeroc,rhsout,one)
    output_spline = zeroc
    call zgemv(charn,nrad,nrad,alph,Aout,nrad,rhsout,one,zeroc,output_spline,one)
    do i=1,nrad
     w_st(i,s) = cmplx(dbsval(rsmall(i),kord,rknot,nrad,real(output_spline)), &
                        dbsval(rsmall(i),kord,rknot,nrad,aimag(output_spline)))
    enddo
    
    
!    print *,'Finished s', s
!  enddo
 enddo
 print *,'Done inversion'
 deallocate(Amat,rhs,rhsout,Aout)


 ! do s=smin,smax
  !print *,maxval(abs(symbols(s)%submat(:,:))),s,wigner(120)%sub3j(s,-1,-119)
 ! enddo
!  bj = 2.0d0
!  ell = 120
  !print *,nr,ns
!  allocate(polR(nr,ns), polH(nr,ns), polT(nr,ns))


!  print *,'here'
! do i=-120,120
!  print *,Zmatrix(4,4)%sub!(-1,i)
! enddo
! print *,'finished'
  !call construct_frequencies(120)
  !do i=-120,120
  !print *, i,nus_mn(i,0), nus_mn(i,1), nus_mn(i,2)
  !enddo

 
! stop
!  call compute_kernels(aj,bj,aj,bj, polR, polH, polT)
  !print *,maxval(abs(polR)),minval(abs(polR))
  !call writefits('radius.fits',r,nr,1,1)
  !call writefits('toroidal_ell_120.fits',polT,nr,ns,1)

 !print *,ned(120.d0,20.d0,120.d0,119.d0,1.d0,120.d0)
 !temv = 0.d0
 !call DRC3JJ (120.d0, 120.d0, -101.d0, 117.d0, lmi, lma, temv, 80, IER)
! print *, temv, ier, lmi, lma
! print *,lmi, lma
! do i=1,80
! print *,i,temv(i)
! enddo
! stop
!  do s = smin, 10!smax
!  temv = 0.d0
!ell = 15
!  do t = -smax, smax 
!   do em = -ell, ell
!     call drc3jj(dble(ell),dble(ell),dble(-em-t),dble(em),lmi,lma,temv,smax,IER) 
!     print *,em,t,lmi,lma
!     if (em == -2 .and. t == 3) then
!       do i=1,smax
!          print *,i,temv(i), wig3j(dble(ell),dble(lmi+i-1),dble(ell),dble(em),dble(t),dble(-em-t))
!       enddo
!       stop
!     endif
!     wigner3j(int(lmi):int(lma),em,t) = temv(1:int(lma-lmi+1))
     !print *,wigner3j(t,em,s),t,em,s
!    enddo
!   enddo
!  enddo
  !print *,minval(abs(wigner3j))
! stop
!  call writefits('Wigner_ell_120.fits',wigner3j,2*ell,2*ell,ns)


!  print *,nr
!   aj = 4.0d0
!   bj = 2.0d0
!   cj = 5.0d0
!   am = 3.0d0
!   bm = -1.0d0
!   cm = -4.0d0
!   cg = 0.0

!   cg = wig3j(aj,bj,cj,am,bm,cm)
!   cg = bcoef(1.0d0,4.d0,5.d0,2.d0,-1.0d0)
   !cg = wig3j(aj,bj,cj,am,bm,cm)
!   print *,cg,bcoef(0.0d0,3.d0,2.d0,2.d0,1.0d0)
!   print *,cg,factorial(5.0d0)!,aj,bj,cj,am,bm,cm
!   call ned(bj,cj,aj,bm,-cm,-am,cg)
   !print *,cg
!   call ned(cj,aj,bj,-cm,am,-bm,cg)
   !print *,cg

!   call wig3j(aj,bj,cj,am,bm,cm,cg)
   !print *,cg
! Verified for a few integer and non-integer values


end program inversion
