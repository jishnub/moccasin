PROGRAM SETUP

 use f90getopt
 implicit none
 include 'params.i'

 integer ell, year, t, s, ellp,int_nyears, int_yearlength, int_yearend, int_yearnum
 integer elmin, elmax, yst,esst, j, ell2, dell
 real*8 nyears, yearlength, yearend, yearnum

 character*3 ellc, ellc2
 character*2 ynum, tc, sc
 character*80 outdir

 logical data_process

  type(option_s):: opts(9)
   opts(1) = option_s( "compute_norms", .false., 'a' )
   opts(2) = option_s( "dell",  .true.,  'b' )
   opts(3) = option_s( "yearnum", .true.,  'c')
   opts(4) = option_s( "lmin", .true.,  'd')
   opts(5) = option_s( "nyears", .true.,  'e')
   opts(6) = option_s( "help", .false.,  'f')
   opts(7) = option_s( "instrument", .true.,  'g')
   opts(8) = option_s( "lmax", .true.,  'h')
   opts(9) = option_s( "kernels", .false., 'i')

   compute_norms = .false.
   yearnum = -1
   elmin = -1
   dell = -1
   yearlength = -1
   elmax = -1
   instrument = 'FFF'
   data_process = .true.

    do
        select case( getopt( "abcdefhi:", opts ) ) ! opts is optional (for longopts only)
            case( char(0) )
                exit
            case( 'a' )
                compute_norms = .true.
                dell = 0
            case( 'b' )
                read(optarg,*) dell 
            case( 'c' )
                read(optarg,*) yearnum
            case( 'd' )
                read(optarg,*) elmin
            case( 'e' )
                read(optarg,*) yearlength
            case( 'f' )
                print *,'Need to provide options of the form'
             print*,' ./setup --lmin 50 --lmax 60 --nyears 5 --yearnum 2 --instrument HMI --dell 2'
                print *,'Note that dell has to be greater than or equal to 0'
                print *,'Optional inputs are --compute_norms (default is not to compute normalization), which forces dell = 0'
                print *,'and --kernels to compute leakage kernels (default is for data processing)'
                print *,'Quitting.'
                stop
            case( 'g' )
                if (optarg == 'HMI' .or. optarg == 'MDI') instrument = optarg
                if (optarg == 'mdi' .or. optarg == 'hmi') instrument = upper(optarg)
            case( 'h' )
                read(optarg,*) elmax
            case( 'i' )
                data_process = .false.

        end select
    end do

    if (compute_norms) dell = 0
    if (yearlength < 0 .or. elmin < 0 .or. dell < 0 .or. yearnum < 0 &
        .or. instrument == 'FFF' .or. elmax < 0) then
     print *,'Need lmin, lmax, dell, nyears, start year, maxell and instrument. Quitting.'
     stop
    endif
 yearend = yearnum + yearlength - 1 

 if (compute_norms) dell = 0
 if (instrument == 'HMI') outdir = outhmidir
 if (instrument == 'MDI') outdir = outmdidir
 call system('mkdir -p '//trim(adjustl(outdir)))
 if (data_process) then
  call system('mkdir -p '//trim(adjustl(outdir))//'/processed')

  int_yearnum = nint((yearnum-1.0)*5.0) + 1
  int_yearlength = nint(yearlength*5.0)
  int_yearend = nint((yearend - 1.0)*5.0)+1

  do year = int_yearnum,int_yearend,int_yearlength
   write(ynum,'(I2.2)') year!nint((year-1)*5.0) + 1
  !print *,year,ynum
 !stop
 !cycle
   do ell = elmin, elmax
!    do ell2 = ell+1,lmax
     ell2 = ell + dell
     if (ell2 > elmax) cycle
     write(ellc,'(I3.3)') ell
     write(ellc2,'(I3.3)') ell2
     
     call system('touch '//trim(adjustl(outdir))//'/processed/'//ellc//'_'//ellc2//'_'//ynum)

   enddo
  enddo
 else 
  call system('mkdir -p '//trim(adjustl(outdir))//'/kernprocess')
  do t = 1, smax
   write(tc,'(I2.2)') t
   esst=max(floor(t/2.d0)*2+1,5)
   do s = esst, smax, 2
    write(sc,'(I2.2)') s
    call system('touch '//trim(adjustl(outdir))//'/kernprocess/'//tc//'_'//sc)
   enddo
  enddo
 endif

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


END PROGRAM SETUP
