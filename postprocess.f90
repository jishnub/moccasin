PROGRAM POSTPROCESS


 implicit none

 integer ell, n, year, yearlength, nsig, nyears, ind, s, t, lmin, lmax, indnl, k
 integer nord, ords(0:25), i, ns, freqmax, ierr, dst, smax, smin, nmodes, freqmin, indal

 real*8 dnu, useless, rsun, cons, speednorm, noisecon
 parameter(rsun = 695.98e6, speednorm = 43692898.859838411)
 character*80 filename, metadata, noisename, varname
 character*3 lch, instrument, track
 character*2 ynum,ynum2, dep
 character*1 ordc, ordcm

 logical exists, sound

 integer, allocatable, dimension(:,:) :: inds
 real*8, allocatable, dimension(:,:) :: alphanl, alphas, temp4, noisecoef, varst, var
 real*8, allocatable, dimension(:,:,:) :: temp2, temp3, noisest, vartemp
 real*8, allocatable, dimension(:,:,:,:) :: temp,noisefreq, noisestfreq, noise
 complex*16, allocatable, dimension(:,:,:,:) :: pow, powst
 !complex*16, allocatable, dimension(:,:,:,:,:) ::

 complex*16 eye
 parameter(eye = (0.d0, 1.d0))
! real*8 ff(0:30,-30:30)


 call getarg(1, instrument)
 call getarg(2, dep)
 call getarg(3, track)
 call getarg(4, ordc)

 sound = .false.
 if (ordc == 't') sound = .true.

 !instrument = 'HMI'


 smax = 50
 smin = 1
 if (.not. sound) smin = 5
 ns = (smax+1)**2

 yearlength = 1

 dnu = 1e6/(360.d0*24.d0*3600.d0 * yearlength)
 lmin = 10
 lmax = 191
 nyears = 14
 if (instrument == 'HMI') then 
  lmin = 10
  nyears = 3
 endif
 freqmin = 1!max(floor(0./dnu),1)
 freqmax = 95!min(floor(2.5/dnu),75)
 nsig = 95!75
 
 allocate(temp(1:ns,0:30,nsig,2),pow(1:freqmax,0:smax,-smax:smax,nyears),&
     powst(1:freqmax,0:smax,-smax:smax,nyears),inds(0:lmax,0:30),&
     noise(1:freqmax,0:smax,-smax:smax,nyears),noisest(0:smax,-smax:smax,nyears),&
     !noisefreq(1:freqmax,0:smax,-smax:smax,nyears),&
     noisestfreq(1:freqmax,0:smax,-smax:smax,nyears),&
     var(0:smax,-smax:smax), varst(0:smax,-smax:smax))!,&
!     temp2(1:ns,0:25,nsig))!, alphanl(lmin:lmax,0:26))

 inds = 0
 nmodes = 0
 if (.not. sound) open(55,file='alpha_'//instrument//'.txt',action='read',position='rewind')
 if (sound) open(55,file='alpha_sound_'//instrument//'.txt',action='read',position='rewind')
 do 
  read(55,*,IOSTAT=ierr) ell,n
  if (ierr .ne.0) exit
  !alphanl(int(ell),int(n)) = useless
  nmodes = nmodes+1
  inds(ell,n) = nmodes
 enddo
 close(55)
 allocate(alphanl(ns,nmodes),noisecoef(ns,nmodes),temp3(nmodes,nsig,2),temp4(nmodes,nsig))

 if (.not. sound) call readfits_normal('alphas_'//instrument//'_tracking_'// &
  track//'_depth_'//dep//'.fits',alphanl,ns,nmodes,1)

 if (sound) call readfits_normal('alphas_sound_'//instrument//'_tracking_'// &
  track//'_depth_'//dep//'.fits',alphanl,ns,nmodes,1)

 do year = 2,4!1,nyears,yearlength
  print *,year
  write(ynum,'(I2.2)') 5*(year-1)+1
  write(ynum2,'(I2.2)') yearlength*5

  do ell = lmin, lmax
   print *,ell
   if (sum(inds(ell,:)) == 0) cycle 
   write(lch,'(I3.3)') ell

   nord = 0
   filename = '/scratch/shravan/'//instrument//'2/'//'tracking'//track//'/bcoef_l_'//lch//'_lp_'//lch//&
   '_year_'//ynum//'_'//ynum2//'.fits'

   metadata = '/scratch/shravan/'//instrument//'2/'//'tracking'//track//'/bcoef_metadata_l_'//lch//'_lp_'//lch//&
   '_year_'//ynum//'_'//ynum2

   noisename = '/scratch/shravan/'//instrument//'2/'//'tracking'//track//'/noise_l_'//lch//'_lp_'//lch//&
   '_year_'//ynum//'_'//ynum2//'.fits'

   varname = '/scratch/shravan/'//instrument//'2/'//'tracking'//track//'/noisevar_l_'//lch//'_lp_'//lch//&
   '_year_'//ynum//'_'//ynum2//'.fits'



   open(1555,file=trim(adjustl(metadata)),status='old',action='read')
   do  
    read(1555,*,IOSTAT=ierr) i,useless,useless
    if (ierr .ne.0) exit
    nord = nord + 1
    ords(nord) = i
   enddo
   close(1555)

   inquire(file=trim(adjustl(filename)), exist=exists)

   if (exists) then

    if (allocated(temp)) deallocate(temp)
    allocate(temp(1:ns,1:nord,1:nsig,2))

    if (allocated(temp2)) deallocate(temp2)
    allocate(temp2(1:ns,1:nord,1:nsig))

    if (year == 1) then
     if (allocated(vartemp)) deallocate(vartemp)
     allocate(vartemp(1:ns,1:nord,1:nsig))
    endif
     

    call readfits(trim(adjustl(filename)),temp(1,1,1,1),ns,nord,nsig)
    call readfits_normal(trim(adjustl(noisename)),temp2(1,1,1),ns,nord,nsig)
    if (year==1) call readfits_normal(trim(adjustl(varname)),vartemp(1,1,1),ns,nord,nsig)

    if (sum(abs(temp)**2) > 1.0) cycle
    do i = 1, nord
     if (ords(i) .lt. 0) cycle
     indnl = inds(ell, ords(i)) 
      
     if (indnl ==0) cycle
     !if (alphanl(ell,ords(i))==0) cycle

     do t = -smax, smax
      if (t==0) cycle
      do s = max(abs(t),smin), smax
       k = mod(s,2)
       if (k == 0 .and. (.not. sound)) cycle
       if (k == 1 .and. sound) cycle

       dst = s-abs(t)
       n = smax - dst - max(smin - dst -1,0)
       ind = s*(s+1) + t + 1
       indal = ind
!       if (k == 0 .and. dst >0) indal = (s-1)*s + t + 1
!       if (k == 0 .and. dst ==0) indal = (s-1)*s + t 

!       if (alphanl(ind,indnl) == 0) cycle
       !if (k == 0) ind = s*(s-1) + t + 1
       cons =  rsun * alphanl(indal,indnl) *(s*(s+1))**0.5
       noisecon =  cons 

       pow(freqmin:freqmax,s,t,year) = pow(freqmin:freqmax,s,t,year) + (temp(ind,i,freqmin:freqmax,1) + &
          eye * temp(ind,i,freqmin:freqmax,2)) * cons 


     noise(freqmin:freqmax,s,t,year) = noise(freqmin:freqmax,s,t,year) + temp2(ind,i,freqmin:freqmax) * cons**2
     if (year==1) var(s,t) = var(s,t) + sum(vartemp(ind,i,freqmin:freqmax)) * cons**4

! THIS WAS THE ORIGINAL SUM   
       powst(freqmin:freqmax,dst,t,year) = powst(freqmin:freqmax,dst,t,year) + (temp(ind,i,freqmin:freqmax,1) + &
          eye * temp(ind,i,freqmin:freqmax,2)) * cons 


       noisestfreq(freqmin:freqmax,dst,t,year) = noisestfreq(freqmin:freqmax,dst,t,year) + temp2(ind,i,freqmin:freqmax) * cons**2
       noisest(dst,t,year) = noisest(dst,t,year) + sum(temp2(ind,i,freqmin:freqmax)) * cons**2
       if (year==1) varst(dst,t) = varst(dst,t) + sum(vartemp(ind,i,freqmin:freqmax)) * cons**4

!       noisefreqst(freqmin:freqmax,dst,t,year) = noisefreqst(freqmin:freqmax,dst,t,year) + (temp(ind,i,freqmin:freqmax,1) + &
!          eye * temp(ind,i,freqmin:freqmax,2)) * noisecon!s*(1.-2*modulo(ell,2)) !* (s * (s+1))**0.5 !ell,ords(i))

      enddo
     enddo
    enddo
   endif

  enddo
 enddo

 
 if (.not. sound) then
  call writefits('tracking'//track//'/power_'//instrument//'_'//dep//'.fits',sum(abs(pow(freqmin:freqmax,:,:,:))**2,1), &
                smax+1,2*smax+1,nyears,1)

  call writefits('tracking'//track//'/power_freq_'//instrument//'_'//dep//'.fits',&
              sum(abs(pow(freqmin:freqmax,:,:,:))**2,4)/nyears, &
                freqmax-freqmin+1,smax+1,2*smax+1,1)

 call writefits('tracking'//track//'/noise_freq_'//instrument//'_'//dep//'.fits',&
                sum(noise(freqmin:freqmax,:,:,:),4)/nyears, &
                freqmax-freqmin+1,smax+1,2*smax+1,1)

  call writefits('tracking'//track//'/noise_'//instrument//'_'//dep//'.fits',&
                sum(noise(freqmin:freqmax,:,:,:),1), &
                smax+1,2*smax+1,nyears,1)

  call writefits('tracking'//track//'/var_'//instrument//'_'//dep//'.fits',var/nyears, &
                smax+1,2*smax+1,1,1)

  call writefits('tracking'//track//'/power_st_'//instrument//'_'//dep//'.fits',sum(abs(powst(freqmin:freqmax,:,:,:))**2,1), &
                smax+1,2*smax+1,nyears,1)


  call writefits('tracking'//track//'/power_st_freq_'//instrument//'_'//dep//'.fits',sum(abs(powst)**2,4)/nyears, &
                freqmax-freqmin+1,smax+1,2*smax+1,1)

  call writefits('tracking'//track//'/noise_st_freq_'//instrument//'_'//dep//'.fits',sum(noisestfreq,4)/nyears, &
                freqmax-freqmin+1,smax+1,2*smax+1,1)


! call writefits('power_st_'//instrument//'_'//dep//'.fits',sum(real(powst(freqmin:freqmax,:,:,:)),4)/nyears, &
!                freqmax-freqmin+1,smax+1,2*smax+1,1)

! call writefits('noise_st.fits',sum(abs(noisefreqst(freqmin:freqmax,:,:,:))**2,1), &
!                smax+1,2*smax+1,nyears,1)

  call writefits('tracking'//track//'/noise_st_'//instrument//'_'//dep//'.fits',noisest, &
                smax+1,2*smax+1,nyears,1)

  call writefits('tracking'//track//'/var_st_'//instrument//'_'//dep//'.fits',varst/nyears, &
                smax+1,2*smax+1,1,1)

 
  else
 
  call writefits('tracking'//track//'/power_sound_'//instrument//'_'//dep//'.fits',sum(abs(pow(freqmin:freqmax,:,:,:))**2,1), &
                smax+1,2*smax+1,nyears,1)

  call writefits('tracking'//track//'/power_sound_freq_'//instrument//'_'//dep//'.fits',sum(abs(pow(:,:,:,:))**2,4)/nyears, &
                freqmax-freqmin+1,smax+1,2*smax+1,1)

! call writefits('noise.fits',sum(abs(noisefreq(freqmin:freqmax,:,:,:))**2,1), &
!                smax+1,2*smax+1,nyears,1)
  !call writefits('tracking'//track//'/noise_sound_'//instrument//'_'//dep//'.fits',noise(:,:,1), &
   !             smax+1,2*smax+1,1,1)

  call writefits('tracking'//track//'/power_sound_st_'//instrument//'_'//dep//'.fits',sum(abs(powst(freqmin:freqmax,:,:,:))**2,1), &
                smax+1,2*smax+1,nyears,1)


 !THIS WAS THE ORIGINAL OUTPUT
  call writefits('tracking'//track//'/power_sound_st_freq_'//instrument//'_'//dep//'.fits',sum(abs(powst(:,:,:,:))**2,4)/nyears, &
                freqmax-freqmin+1,smax+1,2*smax+1,1)


! call writefits('power_st_'//instrument//'_'//dep//'.fits',sum(real(powst(freqmin:freqmax,:,:,:)),4)/nyears, &
!                freqmax-freqmin+1,smax+1,2*smax+1,1)

! call writefits('noise_st.fits',sum(abs(noisefreqst(freqmin:freqmax,:,:,:))**2,1), &
!                smax+1,2*smax+1,nyears,1)

  call writefits('tracking'//track//'/noise_sound_st_'//instrument//'_'//dep//'.fits',noisest(:,:,1), &
                smax+1,2*smax+1,1,1)

  endif
 deallocate(alphanl, inds, pow, powst, temp, temp2, noise, noisest, var, varst)
Contains

!================================================================================

	 SUBROUTINE readfits(filename,readarr,dim1,dime2,dim3)

	  implicit none
	  integer status,unit,readwrite,blocksize,naxes(4)
	  integer group,firstpix, dim3 ,dime2,dim1
	  integer nelements
	  real*8 readarr(dim1,dime2,dim3,2)
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
           naxes(4) = 2
           nelements=naxes(1)*naxes(2)*naxes(3)*2
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

	 SUBROUTINE readfits_normal(filename,readarr,dim1,dime2,dim3)

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
	   
	   !print *,minval(readarr),maxval(readarr)
	  call ftclos(unit, status)
	  call ftfiou(unit, status)


	  end SUBROUTINE readfits_normal


!================================================================================


        SUBROUTINE writefits(filename,temp,nx, ny, nz,nd)

        implicit none 
        integer blocksize,bitpix,naxes(4),unit1,nz , nx, ny,nd
        integer status1,group,fpixel,flag, nelements
	  real*8 temp(nx,ny,nz,nd)
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
	  naxes(4)=nd
	  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
	  extend=.false.
	  group=1
	  fpixel=1

	  call ftphpr(unit1,simple,bitpix,4,naxes,0,1,extend,status1)
	  call ftpprd(unit1,group,fpixel,nelements,temp,status1)
	  call ftclos(unit1, status1)
	  call ftfiou(unit1, status1)


	 end SUBROUTINE writefits
!================================================================================


END PROGRAM POSTPROCESS
