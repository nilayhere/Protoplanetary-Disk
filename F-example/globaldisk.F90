!=========================================================================================
!
!       MAIN PROGRAM TO EXECUTE FULL GLOBAL DISK DUST EVOLUTION
!
!=========================================================================================

#include "global.h"
#include "parameter.h"

program GlobalDisk

   use, intrinsic :: iso_fortran_env
   implicit none

   include 'mpif.h'
   
   integer :: numtasks, rank, ierr, rc, len, e
   character*(MPI_MAX_PROCESSOR_NAME) name

   type :: zdata
      
      double precision :: n
      double precision :: rho_g, temperature, alpha, c_s, zvalue, rho_d
      double precision :: hist(NH), delz

   end type zdata

   type(zdata), dimension(:), allocatable, target :: z

   integer :: i, j, row, col, nz, un, nst, ii, iter, ttt
   integer :: print_counter, time_counter, step_flag

   double precision :: rhog_0, collmass, radius, H, Omega, n_new
   double precision :: totaltime, del_time
   double precision :: hist_bin(NH+1), bin_mass(NH)
   double precision :: a(NH), a_edge(NH+1), print_time, cd_time
   double precision :: array3(NH), mass_frac(NH), num_frac(NH)

   double precision, dimension(:,:), allocatable :: param
   double precision, dimension(:,:), allocatable :: print_before_s, print_after_s
   double precision, dimension(:,:), allocatable :: mass_array_in, mass_array_out

   double precision, dimension(:), allocatable :: z_grid, cs_array, alpha_array
   double precision, dimension(:), allocatable :: gas_den_array, z_array

   character (len=100) :: cwd
   character (len=200) :: filename
   character (len=20)  :: ft, ft1, ft2, X1, X2

!====================================================================================
!
!                       MPI BLOCK : INITIALIZE & SET PROCESSORS
!
!====================================================================================

!  Initialize MPI library

   call MPI_INIT(ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print*, 'ERROR STARTING MPI LIBRARY : TREMINATING '
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if

!  Get the number of processors the job will be using

   call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

!  Get the rank of the processorthe thread is running on :
!  Every processor is assigned a different rank starting with 0

   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

!  Get the name of the processor (Usually the hostname)

   call MPI_GET_PROCESSOR_NAME(name, len, ierr)

   if (ierr .ne. MPI_SUCCESS) then
      print*, 'ERROR GETTING PROCESSOR NAME : TERMINATING JOBS'
      call MPI_ABORT(MPI_COMM_WORLD, rank, ierr)
   end if

!  print*, 'processor Name : ', rank

!  MPI LIBRARY IS INITIALIZED AND SETUP IS COMPLETE. 

!====================================================================================
!
!                     INPUT GRID : 
!
!==================================================================================== 

! Get the current working directory: cwd

   call GETCWD(cwd)

! Read single column parameters:

   allocate(param(NR,6))

   un=22+rank
   open(unit=un, file='parameter.inp', status='old', action='read')
   
   do row=1,NR
      read(un,*) (param(row,col), col=1,6)
   end do

   close(un)

   nz=INT(param(rank+1,1))
   rhog_0=param(rank+1,2)
   collmass=param(rank+1,3)
   radius=param(rank+1,4)*AU
   H=param(rank+1,5)
   Omega=param(rank+1,6)

   deallocate(param)

   allocate(z(nz))

   print*, rank,nz, rhog_0, collmass, radius, H

! Check for interim files first: 

   allocate(param(nz,NH+7))

   ft='(I2.2)'
   write(X1,ft) rank+1

   print*, X1,rank+1

   filename = trim(adjustl(cwd))//'/interim/interim_'//trim(adjustl(X1))//'.dat'

   if (access(filename,' ') .eq. 0) then
      un=111+rank
      open(unit=un, file=filename, status='old', action='read')
      
      do row=1,nz
         read(un,*) (param(row,col), col=1,NH+7)
      end do

      close(un)

   else

      filename = trim(adjustl(cwd))//'/initfiles/start_'//trim(adjustl(X1))//'.dat'

      un=111+rank
      open(unit=un, file=filename, status='old', action='read')

      do row=1,nz
         read(un,*) (param(row,col), col=1,NH+7)
      end do

      close(un)

   end if
   
   do row=1,nz
      z(row)%hist=param(row,1:NH)
      z(row)%rho_g=param(row,NH+1)
      z(row)%rho_d=param(row,NH+2)
      z(row)%n=param(row,NH+3)
      z(row)%alpha=param(row,NH+4)
      z(row)%temperature=param(row,NH+5)
      z(row)%zvalue=param(row,NH+6)
      z(row)%delz=param(row,NH+7)
   end do
   
   deallocate(param)

!  Read edge values for vertical cells 

   allocate(z_grid(nz+1))

   allocate(param(NR,NZMAX+1))

   un=204+rank
   open(unit=un, file='zgrid.inp', status='old', action='read')
   do row=1,NR
      read(un,*) (param(row,col), col=1,NZMAX+1)
   end do
   close(un)

   z_grid=param(rank+1,1:nz+1)

   deallocate(param)

!  Read histogram bin: equispaced in log scale

   allocate(param(NH+1,2))

   un=55+rank
   open(unit=un, file='hist_bin.inp', status='old', action='read')
   do row=1,NH+1
      read(un,*) (param(row,col), col=1,2)   
   end do
   hist_bin=param(:,1)
   bin_mass=param(1:NH,2)
   deallocate(param)

!  Read dust size

   allocate(param(NH+1,2))
   un=306+rank
   open(unit=un, file='dust_size.inp', status='old', action='read')
   do row=1,NH+1
      read(un,*) (param(row,col), col=1,2)
   end do
   a=param(1:NH,2)
   a_edge=param(:,1)
   deallocate(param)

!  Construct array for sound speed, gas density, alpha, z value

   allocate(cs_array(nz))
   allocate(alpha_array(nz))
   allocate(gas_den_array(nz))
   allocate(z_array(nz))

   do i=1,nz
      z(i)%c_s=sqrt((KB*z(i)%temperature)/MUMP)
      cs_array(i)=z(i)%c_s
      alpha_array(i)=z(i)%alpha
      gas_den_array(i)=rhog_0*exp(-(z(i)%zvalue**2d0)/(2d0*H*H))
      z_array(i)=z(i)%zvalue
   end do

   
   allocate(mass_array_in(nz,NH))

   allocate(mass_array_out(nz,NH))

!  MAIN CODE STARTS HERE: CALLING COLLISION ROUTINE FOR EACH GRID

   print_counter=1

   totaltime=0d0

!  un=299

!  Open a file 'timeline.dat': This file will contain all the times settling 
!  has been performed.
!  The index number will match the index of the file number.
!  The code will keep on appending the time to this file.

!  open(unit=un, file='timeline.dat', status='unknown', action='write')
!  write(un,*) totaltime/SECYR
!  close(un)

   del_time=0.1d0*SECYR

   step_flag=0

   do ttt=1,1

      do i=1,nz

!        print*, 'STARTING COLLISION ALGORITHM FOR ZONE : ', i,'   ',rank

         if (step_flag .eq. 0) then

            cd_time = 100000d0*SECYR

            nst = 5000

         else

            nst = 1

         end if

!        if (step_flag .gt. 0) then

!           if (rank .lt. 5) then
         
!              nst = 1500

!           else

!              nst = 3000

!           end if

!        end if

!        if (step_flag .le. 2) then
!           cd_time = 5000d0*SECYR
!        end if

         if (step_flag .eq. 3) then
            cd_time = 300d0*SECYR
         end if

         if (step_flag .eq. 4) then
            cd_time = 500d0*SECYR
         end if
         
         if ((step_flag .gt. 4) .and. (step_flag .le. 7))then
            cd_time = 1000d0*SECYR
         end if 
         
         if (step_flag .gt. 7) then
            cd_time = 5000d0*SECYR
         end if

         array3=0d0

         if (z(i)%rho_d .gt. 1d-50) then

         call silicate_collision_SF(rank+1,z(i)%hist, z(i)%rho_d, z(i)%n, z(i)%rho_g, &
                                    z(i)%temperature, z(i)%alpha, &
                                    radius, z(i)%zvalue, H, cs_array(i), cd_time, &
                                    Omega, bin_mass, hist_bin, &
                                    nst, array3, n_new)

         z(i)%n=n_new

         z(i)%hist=array3/sum(array3)

         mass_array_in(i,:)=array3

         else

         mass_array_in(i,:)=z(i)%hist*z(i)%rho_d

         end if

!        print*, 'THE NUMBER OF ZONE COMPLETED FOR COLLISION : ',i,'   ',rank

      end do

!     step_flag = step_flag + 1

      totaltime=totaltime+del_time

      allocate(print_before_s(nz, NH+7))

      do j=1,nz
   
         print_before_s(j,1:NH)=z(j)%hist
         print_before_s(j,NH+1)=z(j)%rho_g
         print_before_s(j,NH+2)=z(j)%rho_d
         print_before_s(j,NH+3)=z(j)%n
         print_before_s(j,NH+4)=z(j)%alpha
         print_before_s(j,NH+5)=z(j)%temperature
         print_before_s(j,NH+6)=z(j)%zvalue
         print_before_s(j,NH+7)=z(j)%delz

      end do


      un=333+rank

      open(unit=un, file=trim(adjustl(cwd))//'/interim/'//'interim_'//trim(adjustl(X1))// &
                     '_'//trim(adjustl(X2))//'.dat', status='unknown', action='write')

      do j=1,nz

         write(un,*) print_before_s(j,:)

      end do

      close(un)


      un=604

      ft1='(I3.3)'

      write(X2,ft1) print_counter

      open(unit=un, file=trim(adjustl(cwd))//'/out/'//'sb_'//trim(adjustl(X1))// &
                        '_'//trim(adjustl(X2))//'.dat', status='unknown', action='write')

      do j=1,nz

         write(un,*) print_before_s(j,:)

      end do

      close(un)

      deallocate(print_before_s)

      if (step_flag .eq. 0) then
         iter = 1000
      elseif (step_flag .eq. 1) then
         iter = 1000
      elseif (step_flag .eq. 2) then
         iter = 3000
      elseif (step_flag .eq. 3) then
         iter = 5000
      elseif ((step_flag .ge. 4) .and. (step_flag .le. 7)) then
         iter = 10000
      elseif (step_flag .gt. 7) then
         iter = 50000
      end if

!     print*, 'STARTING SETTLING AND DIFFUSION ALGORITHM :'

      do j=1,iter !INT(DELTATIME*SECYR/del_time)

         mass_array_out = 0d0

         call settle_the_dust(nz, mass_array_in, a, z_grid, z_array, cs_array, H, rhog_0, &
                           alpha_array, Omega, del_time, mass_array_out)

         mass_array_in = mass_array_out

      end do

!     print*, 'SETTLING AND DIFFUSION ALGORITHM EXECUTED SUCCESSFULLY '

      do j=1,nz
      
         z(j)%hist=mass_array_out(j,:)/sum(mass_array_out(j,:))

         z(j)%rho_d = sum(mass_array_out(j,:))

         if (sum(mass_array_out(j,:)) .ne. 0d0) then

            mass_frac=mass_array_out(j,:)/sum(mass_array_out(j,:))

         else

            mass_frac = 0d0

         end if

         num_frac = 0d0

         do ii=1,NH

            num_frac(ii) = mass_array_out(j,ii)/(10d0**bin_mass(ii))

         end do

         num_frac=num_frac/sum(num_frac)

         if (sum(num_frac) .gt. 0d0) then

            z(j)%n=sum(mass_array_out(j,:))/dot_product(10d0**bin_mass, num_frac)

         else

            z(j)%n = 0d0

         end if

      end do
     
      allocate(print_after_s(nz, NH+7))

      do j=1,nz
   
         print_after_s(j,1:NH)=z(j)%hist
         print_after_s(j,NH+1)=z(j)%rho_g
         print_after_s(j,NH+2)=z(j)%rho_d
         print_after_s(j,NH+3)=z(j)%n
         print_after_s(j,NH+4)=z(j)%alpha
         print_after_s(j,NH+5)=z(j)%temperature
         print_after_s(j,NH+6)=z(j)%zvalue
         print_after_s(j,NH+7)=z(j)%delz

      end do

!  Save final file after settling: AfeterSettling********.dat

      un=302+rank

!     ft2='(I3.3)'
!     write(X2, ft2) print_counter

      open(unit=un, file=trim(adjustl(cwd))//'/out/'//'op_'//trim(adjustl(X1))// &
                     '_'//trim(adjustl(X2))//'.dat', status='unknown', action='write')

      do j=1,nz

         write(un,*) print_after_s(j,:)

      end do

      close(un)

!     un=333+rank

!     open(unit=un, file=trim(adjustl(cwd))//'/interim/'//'interim_'//trim(adjustl(X1))// &
!                    '_'//trim(adjustl(X2))//'.dat', status='unknown', action='write')

!     do j=1,nz

!        write(un,*) print_after_s(j,:)

!     end do

!     close(un)

      deallocate(print_after_s)

      
!  Save total time to the file timeline.dat

!     un=299
      
!     open(unit=un, file='timeline.dat', access='append', status='old')
!     write(un,*) totaltime/SECYR
!     close(un)

      print_counter=print_counter+1

      step_flag = step_flag+1

   end do

   call MPI_FINALIZE(ierr)

end program GlobalDisk


























