
 PROGRAM ANALYSIS

  use f90getopt
  use data_analysis

  implicit none
  integer ell, ellp,i, nmodels
  real*8 tfin, tstart, nyears, yearnum
  character*3 yearc
  character*3 ellc
  character*80 filename

  type(option_s):: opts(9)
   opts(1) = option_s( "compute_norms", .false., 'a' )
   opts(2) = option_s( "flow_analysis",  .true.,  'b' )
   opts(3) = option_s( "yearnum", .true.,  'c')
   opts(4) = option_s( "ell", .true.,  'd')
   opts(5) = option_s( "ellp", .true.,  'e')
   opts(6) = option_s( "nyears", .true.,  'f')
   opts(7) = option_s( "help", .false.,  'h')
   opts(8) = option_s( "instrument", .true.,  'i')
   opts(9) = option_s( "artificial_data", .true.,  'j')

    compute_norms = .false.
    flow_analysis = .true.
    yearnum = -1
    ell = -1
    ellp = -1
    nyears = -1
    instrument = 'FFF'
    nmodels = 1

    do
        select case( getopt( "abcdefhij:", opts ) ) ! opts is optional (for longopts only)
            case( char(0) )
                exit
            case( 'a' )
                print *, 'option compute_norms'
                compute_norms = .true.
            case( 'b' )
                if (optarg == 'F') flow_analysis = .false.
            case( 'c' )
                read(optarg,*) yearnum
            case( 'd' )
                read(optarg,*) ell
            case( 'e' )
                read(optarg,*) ellp
            case( 'f' )
                read(optarg,*) nyears
            case( 'h' )
                print *,'Need to provide options of the form ./analyze --ell 50 --ellp 60 --nyears 5 --yearnum 2 --instrument HMI'
                print *,'Note that ellp has to be greater than or equal to ell'
                print *,'Optional inputs are --compute_norms (default is not to compute normalization'
                print *,'and --flow_analysis (default is true). If you want sound-speed analysis, run with --flow_analysis F'
                print *,'Quitting.'
                stop
            case( 'i' )
                instrument = optarg
                !if (optarg == 'HMI' .or. optarg == 'MDI') instrument = optarg
                !if (optarg == 'mdi' .or. optarg == 'hmi') instrument = upper(optarg)
            case( 'j' )
                read(optarg,*) nmodels
        end select
    end do
    if (nyears < 0 .or. ell < 0 .or. ellp < 0 .or. yearnum < 0 &
        .or. instrument == 'FFF') then
     print *,'Need ell, ellp, nyears, start year and instrument. Quitting.'
     stop
    elseif (ellp < ell) then
     print *,'Need ellp >= ell. Quitting'
     stop
    endif
 !call getarg(1,yearc)
 !read(yearc,*) yearnum

 !call getarg(2,ellc)
 !read(ellc,*) ell

 !call getarg(3,ellc)
 !read(ellc,*) ellp 

! ellp = ell
 !call getarg(4,yearc)
 !read(yearc,*) nyears
 !call readheader
 !stop
! call glob_dir ()
!  stop
 !MODE_ANALYSIS = .TRUE.
 call cpu_time(tstart)
!  ell = 36
!  ellp = ell
! yearnum =  1
!nyears = 1
 call basic_setup
 call analyzethis(nmodels, nyears, yearnum) 
 call cpu_time(tfin)

! write(112, *)
! write(112, *) 'Time taken'
! write(112, *)
! write(112, *) (tfin - tstart)/60.
! write(112, *)
 close(112)

! write(ellc,'(I3.3)') ell
! write(yearc,'(I3.3)') nint(yearnum*10)

! call system('rm /scr28/shravan/'//instrument//'/processed/'//yearc//'_'//ellc)
 !call system('rm /homea/shravan/QDP/lengs_'//ellc//'_'//ellc//'_'//yearc)

 Contains

 Subroutine glob_dir

  implicit none
  integer ell, ierr, j, daynum
  character*3 ellc
  character*4 dayc
  character*100 dir

  do ell=1,29!lmin,lmax
   write(ellc, '(I3.3)') ell
   call system('rm /scr28/shravan/locations/locs_ell_'//ellc)
   call system('show_info -qP hmi.v_sht_gf_72d[]['//ellc//'] > /scr28/shravan/locations/locs_ell_'//ellc)
   !call system('show_info -qP mdi.vw_V_sht_gf_72d[]['//ellc//'] > /scr28/shravan/locations/locs_ell_'//ellc)
   open(22,file='/scr28/shravan/locations/locs_ell_'//ellc,action='read',position='rewind',status='unknown')
   daynum = 6328
   !daynum = 1216
   do j=1,41
   !do j=1,74
    read(22,'(A)',IOSTAT=ierr) dir
    if (ierr .ne.0 .and. j <41) then
    !if (ierr .ne.0 .and. j <74) then
     print *,'File not found', ell, daynum
     exit
    endif
    write(dayc, '(I4.4)') daynum
    print *,trim(adjustl(dir))//'/data.fits'//' /scr28/shravan/HMI/data/HMI_'//ellc//'_'//dayc//'.fits'
    !print *,trim(adjustl(dir))//'/data.fits'//' /tmp29/shravan/MDI/data/MDI_'//ellc//'_'//dayc//'.fits'
    call system('cp '//trim(adjustl(dir))//'/data.fits'//' /scr28/shravan/HMI/data/HMI_'//ellc//'_'//dayc//'.fits')
    !call system('cp '//trim(adjustl(dir))//'/data.fits'//' /tmp29/shravan/MDI/data/MDI_'//ellc//'_'//dayc//'.fits')
    daynum = daynum + 72
   enddo
  enddo

 end subroutine glob_dir

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

end program ANALYSIS
