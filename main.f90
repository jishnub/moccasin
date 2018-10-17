
 PROGRAM MAIN

  use f90getopt
  use kernels

  implicit none
  integer s,t, cond, yearnum
  logical asymptotics, leake, sound_leaks!, compute_wigner
  character*2 sc, tc
  character*3 trackch
  complex*16, allocatable, dimension(:,:,:) :: Lfun

  type(option_s):: opts(7)
   opts(1) = option_s( "asymptotics", .false., 'a' )
   opts(2) = option_s( "leakage", .false.,  'b')
   opts(3) = option_s( "s", .false.,  'c')
   opts(4) = option_s( "t", .false.,  'd')
   opts(5) = option_s( "help", .false.,  'e')
   opts(6) = option_s( "sound_leaks", .false.,  'f')
   opts(7) = option_s( "instrument", .true.,  'g')
  
   sound_leaks = .false.
   leake = .false.
   s = -2000; t = -2000;
   compute_wigner = .true.
   instrument = 'FFF'

    do
        select case( getopt( "abcdefg:", opts ) ) ! opts is optional (for longopts only)
            case( char(0) )
                exit
            case( 'a' )
                asymptotics = .true.
            case( 'b' )
                leake = .true.
                asymptotics = .false.
                compute_wigner = .false.
            case( 'c' )
                read(optarg,*) s
            case( 'd' )
                read(optarg,*) t
            case( 'e' )
                print *,'Need to provide options of the form ./analyze --ell 50 --ellp 60 --nyears 5 --yearnum 2 --instrument HMI'
                print *,'Note that ellp has to be greater than or equal to ell'
                print *,'Optional inputs are --compute_norms (default is not to compute normalization'
                print *,'and --flow_analysis (default is true). If you want sound-speed analysis, run with --flow_analysis F'
                print *,'Quitting.'
                stop
            case( 'f' )
                sound_leaks = .true.
            case( 'g' )
                if (optarg == 'HMI' .or. optarg == 'MDI') instrument = optarg
                if (optarg == 'mdi' .or. optarg == 'hmi') instrument = upper(optarg)

        end select
    end do

   if (instrument == 'FFF') then
     print *,'Need Instrument. Quitting'
     stop
   endif

   print *,'Standard basic setup'
   call basic_setup!(compute_wigner)


   if (.not. leake) &
    call compute_kernels_all(s, ntot, asymptotics)
    

 if (leake) then
  if((s==-2000) .or. (t == -2000)) then
     print *,'Need s and t when leakage is selected. Quitting'
     stop
  else
   write(sc,'(I2.2)') s 
   write(tc,'(I2.2)') t 
   write(trackch,'(I3.3)') int(1000*trackrate)
   
   allocate(Lfun(lmin+6:lmax-6,0:30,1:63))
   Lfun = 0.d0
   yearnum = 1
   cond = 1
   if (sound_leaks) cond = 0 
   call compute_Q3(4,2,4,3,t,s,Lfun,1.d0,yearnum,cond)
  
   call system('mkdir -p /scratch/shravan'//instrument//'/kernels')
   call system('mkdir -p /scratch/shravan'//instrument//'/kernels/tracking'//trackch)

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

  endif
 endif

 call exit()

 Contains

! ---------------------------------
FUNCTION Upper(s1)  RESULT (s2)
CHARACTER(*)       :: s1
CHARACTER(LEN(s1)) :: s2
CHARACTER          :: ch
INTEGER,PARAMETER  :: DUC = ICHAR('a') - ICHAR('A')
INTEGER            :: i

DO i = 1,LEN(s1)
   ch = s1(i:i)
   IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)-DUC)
   s2(i:i) = ch
END DO
END FUNCTION Upper 

! ---------------------------------

end program MAIN
