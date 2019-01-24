MODULE FREQUENCIES

 implicit none
 include 'params.i'

 type :: subcell
  real*8, allocatable, dimension(:,:) :: mn
 end type 
 
 type(subcell) :: nus(lmin:lmax)

 type :: sub1
  real*8, allocatable, dimension(:) :: ords, ems
 end type 
  
 type(sub1) :: freqnu(lmin:lmax), en(lmin:lmax), a1(lmin:lmax),&
  amps(lmin:lmax), fwhm(lmin:lmax), backg(lmin:lmax)!, leakage(0:304)


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
  logical lexist, mode_analysis
  character*80 workdir
  mode_analysis = .false.

   !print *,ll,llp
!  print *,mode_analysis
!  if (mode_analysis) then

  i = getcwd(workdir)
   write(ellc,'(I3.3)') ll
   write(ellpc,'(I3.3)') llp
   write(yc,'(I2.2)') ynum

   inquire(file=trim(adjustl(workdir))//'/lengs_'//ellc//'_'//ellpc//'_'//yc, exist=lexist)
   if (lexist) call system('rm '//trim(workdir)//'/lengs_'//ellc//'_'//ellpc)

  call system("awk '{ print NF }' "//inpfile//" > "//trim(workdir)//"/lengs_"//ellc//"_"//ellpc//'_'//yc)

  open(52,file=trim(workdir)//'/lengs_'//ellc//'_'//ellpc//'_'//yc,&
           status='old',position='rewind',action='read')
  read(52,*) nrecs
  close(52)


  call system("awk 'END {print NR}' "//inpfile//" > "//trim(workdir)//"/lengs_"//ellc//"_"//ellpc//'_'//yc)

  open(52,file=trim(workdir)//'/lengs_'//ellc//'_'//ellpc//'_'//yc,&
           status='old',position='rewind',action='read')
  read(52,*) npoints
  close(52)

  inquire(file=trim(workdir)//'/lengs_'//ellc//'_'//ellpc//'_'//yc, exist=lexist)
  if (lexist) call system("rm "//trim(workdir)//"/lengs_"//ellc//"_"//ellpc//'_'//yc)

!  else
 !  print *,'READING CANONICAL MODE PARAMETER FILE'
 !  if (instrument == 'MDI') then
 !   npoints = 2139
 !  else
 !   npoints = 2077
 !  endif
   !nrecs = 84
 
  !endif
  
  
  allocate(ells(npoints), ensall(npoints), freqs(npoints),amplitudes(npoints),&
        fwhmall(npoints), backgall(npoints), splits(nrecs,npoints), &
        acoefs(1:nsplits,npoints))

  open(334,file=inpfile,position='rewind',status='old')
  do i=1,npoints
   read(334, *) splits(:,i)
  enddo
  close(334)

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

!  do l = lmin,lmax
  do l = ll, llp

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
    if (allocated(a1(l)%ords)) deallocate(a1(l)%ords)

    minord = int(ensall(ind_data(1)))
    maxord = int(ensall(ind_data(ndata)))

!     print *,minord,maxord,l,ensall(ind_data(1:ndata))
!     print *,minord,maxord,l, nsplits,ensall(ind_data(1)),ensall(ind_data(ndata)),ndata
    allocate(freqnu(l)%ords(minord:maxord), en(l)%ords(minord:maxord), amps(l)%ords(minord:maxord), &
      fwhm(l)%ords(minord:maxord), backg(l)%ords(minord:maxord), a1(l)%ords(minord:maxord),&
      adata(1:nsplits,minord:maxord),Po(-l:l, 1:Nsplits),indicesn(1:ndata))
!Po(-l:l, 1:(Nsplits/2+1))
    
      en(l)%ords = -1.d0
      freqnu(l)%ords = 0.d0
      amps(l)%ords = 0.d0
      fwhm(l)%ords = 0.d0
      backg(l)%ords = 0.d0
      a1(l)%ords = 0.d0
      adata = 0.d0

    indicesn = int(ensall(ind_data(1:ndata)))
    en(l)%ords(indicesn) = ensall(ind_data(1:ndata))
    freqnu(l)%ords(indicesn) = freqs(ind_data(1:ndata))
    amps(l)%ords(indicesn) = amplitudes(ind_data(1:ndata))
    fwhm(l)%ords(indicesn) = fwhmall(ind_data(1:ndata))
    backg(l)%ords(indicesn) = backgall(ind_data(1:ndata))
    !if (lbound(freqnu(l)%ords,1) .ne. minval(en(l)%ords)) &
    ! print *,minord,maxord,lbound(freqnu(l)%ords),minval(en(l)%ords)
    adata(:,indicesn) = acoefs(:,ind_data(1:ndata))
    a1(l)%ords(indicesn) = acoefs(1,ind_data(1:ndata))
    Po(-l:l,:) = Pf(l, Nsplits)

    ordmin = int(en(l)%ords(minord - 1 +minloc(abs(freqmin - freqnu(l)%ords(:)),1)))
    ordmax = int(en(l)%ords(minord-1+minloc(abs(freqmax - freqnu(l)%ords(:)),1)))

    if (allocated(nus(l)%mn)) deallocate(nus(l)%mn)
    allocate(nus(l)%mn(-l:l,minord:maxord))

    do i=minord,maxord
     nus(l)%mn(:,i) = freqnu(l)%ords(i) 
     !do j=1,nsplits/2

     do j=-l,l ! SYNODIC CORRECTION TO THE ROTATION RATE - ONLY THE a1 COEFFICIENT IS ALTERED
       nus(l)%mn(j,i) = nus(l)%mn(j,i) + (adata(1,i) - rotEarth)* Po(j,1) !(trackrate - rotEarth)* Po(j,1)!
     enddo

     do j=2,nsplits ! NOW ADDING THE REST OF THE SPLITTING COEFFICIENTS
      nus(l)%mn(:,i) = nus(l)%mn(:,i) + adata(j,i) * Po(:,j)!a1_data(i)*Po(:,1) + a3_data(i)*Po(:,2) &
             !+ a5_data(i)*Po(:,3)
     enddo

    enddo
    deallocate(adata)

   enddo

 
  deallocate(acoefs, ensall, freqs, amplitudes, fwhmall, backgall, Po)
end subroutine construct_frequencies


!------------------------------------------------------------------------
END MODULE FREQUENCIES
