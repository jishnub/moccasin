
 PROGRAM ANALYSIS

  use data_analysis

  implicit none
  integer ell, ellp, nyears, yearnum
  real*8 tfin, tstart
  character*2 yearc
  character*3 ellc

 call getarg(1,yearc)
 read(yearc,*) yearnum

 call getarg(2,ellc)
 read(ellc,*) ell

 call getarg(3,ellc)
 read(ellc,*) ellp 

 call getarg(4,yearc)
 read(yearc,*) nyears
 !call readheader
 !stop
 MODE_ANALYSIS = .TRUE.
 call cpu_time(tstart)
 call analyzethis(yearnum, ell, ellp, nyears) 
 call cpu_time(tfin)

! write(112, *)
! write(112, *) 'Time taken'
! write(112, *)
! write(112, *) (tfin - tstart)/60.
! write(112, *)
 close(112)

 write(ellc,'(I3.3)') ell
 write(yearc,'(I2.2)') yearnum

 call system('rm /scratch/shravan/'//instrument//'/processed/'//yearc//'_'//ellc)
 call system('rm /home/shravan/QDP/lengs_'//ellc//'_'//ellc//'_'//yearc)

end program ANALYSIS
