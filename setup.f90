PROGRAM SETUP

 implicit none

 integer ell, year, yearlength, t, s, smax, lmin, lmax, smin, yst,esst

 character*2 ynum, tc, sc
 character*3 ellc, instrument

 logical data_process, compute_noise

 data_process = .false.
 instrument = 'HMI'
 compute_noise = .true.
!----------------------------

 yst =2
 yearlength = 14

 smin = 4
 smax = 30
 lmin = 30
 lmax = 249

 if (instrument =='HMI') then 
  lmin = 50
  yearlength = 6
 endif
 !if (compute_noise) then
 ! yst = 1
 ! yearlength = 6
 !endif
 
 if (data_process) then
  do year = yst,yearlength
   write(ynum,'(I2.2)') year
   do ell = lmin, lmax
    write(ellc,'(I3.3)') ell
  !   call system('touch /scratch/shravan/MDI/processed/'//ynum//'_'//ellc)
     !call system('touch /scratch/shravan/HMI/processed/'//ynum//'_'//ellc)
     call system('touch /scratch/shravan/'//instrument//'/processed/'//ynum//'_'//ellc)
   enddo
  enddo
 else 
  do t = 1, smax
   write(tc,'(I2.2)') t
   esst=max(floor(t/2.d0)*2+1,5)
   do s = esst, smax, 2
    write(sc,'(I2.2)') s
    call system('touch /scratch/shravan/kernprocess/'//tc//'_'//sc)
   enddo
  enddo
 endif



END PROGRAM SETUP
