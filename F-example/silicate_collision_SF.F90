!---------------------------------------------------------------------------------------------------------
!
!                                 SUBROUTINE silicate_collision_SF( )
!
!  This subtoutine implements the collision model of Sticking and Frgamentation 
!  The particles are all taken as silicate grains with material density 3.4 g/cm^3
!  
!  The subroutine runs for a preset time interval (delta_time), passed by the main program
!  and returns a final dust distribution as an output.
!  
!  
!  The sticking threshold velocity and fragmentation threshold velocity does not depend on this 
!  routine and is passed by main program. So, any collision model with S and F can be implemented 
!  using this routine. The threshold velocities should be changed in the header file.
!
!---------------------------------------------------------------------------------------------------------
!
!                          : ARGUMENTS
!
!              in_arr: Input file with mass in each bin
!      cell_dust_nden: Total dust number density in the cell
!        cell_gas+den: Gas density of the particular cell
!           cell_temp: Temperature of the cell
!           cel_alpha: Turbulence parameter \alpha of the cell
!            coll_rad: Distance of the column from the central star (cm)
!              cell_z: Height of the cell above midplane (cm)     
!            coll_sch: Scale-height of the column
!             cell_cs: Thermal gas velocity of the cell
!          delta_time: Time span through which the 0^D simulation will run
!                   t: Initial time when the collision model is initiated
!          coll_omega: Angular frequency of the vertical column
!                 bin: array of mass values of histogram bin (NH), log_space
!            bin_edge: Edge values of the bin (NH+1)
!             out_arr: Final dust distribution of the cell (output)
!        dust_den_out: Resulting dust number density (output)
!
!
!--------------------------------------------------------------------------------------------------------



#include "global.h"
#include "parameter.h"


subroutine silicate_collision_SF(rank1,in_arr,  dustmass, cell_dust_nden, cell_gas_den, cell_temp, cell_alpha, coll_rad, &
                                  cell_z, coll_sch, cell_cs, delta_time, coll_omega, bin, &
                                    bin_edge, nstep, out_arr, dust_den_out)

   implicit none

!  LOCAL AND PASSING VARIABLES.

   integer :: ii,j,kk,ll,mm,nn,pp,N,i,nstick,n_coll, nnn,counter,un,tracker, coll_flag

   integer, intent(in) :: nstep, rank1
   double precision, intent(in) :: delta_time, in_arr(NH), dustmass
   double precision, intent(out) :: out_arr(NH), dust_den_out
   double precision, intent(in) :: cell_dust_nden, cell_gas_den, cell_temp, cell_alpha, coll_omega
   double precision, intent(in) :: coll_rad, cell_z, coll_sch, cell_cs, bin_edge(NH+1), bin(NH)
   double precision, dimension(:), allocatable :: c1, c2, wc1, wc2, vc
   double precision, dimension(:), allocatable :: nc1, nc2, wnc1, wnc2, vnc
   double precision, dimension(:), allocatable :: s1, s2, ws1, ws2, vs, ss1, ss2, wss1, wss2
   double precision, dimension(:), allocatable :: ns1, ns2, wns1, wns2, vns, nss1, nss2, wnss1, wnss2
   double precision, dimension(:), allocatable :: f1, f2, wf1, wf2, vf
   double precision, dimension(:), allocatable :: e1, e2, we1, we2, ve
   double precision, dimension(:), allocatable :: tt1, tt2, tt3
   double precision :: bin_num_den(NH), rn1(2), bin_wid_rs(NH)
   integer :: nh1, nh2
   integer :: num_array(NH), frc1(NH), frc2(NH), frc3(NH),frc4(NH), numb1(NH), numb2(NH), numb(NH)
   double precision, dimension(:), allocatable :: hs, hns1, hns2, non_zero_array
   double precision, dimension(:), allocatable :: array1, array2, par_frac_1, par_frac_2, v, randarray
   double precision, dimension(:), allocatable :: weight1, weight2, coll_time, prob1,prob2,prob3,prob4,prob
   double precision :: bin_frac(NH), weights(NH), t_coll_min, t2,t3, pit, t
   double precision :: w_from_f(NH), bin_num(NH), tcollmin, r1, r2
   double precision :: temp_frag_mass(NH), max_array(NH)
!  double precision :: n_from_f(NH), temp_frag_num(NH), temp_erod_num(NH)
   double precision :: dust_num_den, m_in_bin, w_final(NH), crosssection, l2, erod_max
   double precision, dimension(:), allocatable :: tc, vcc, lc, newnden
   double precision, dimension(:,:),allocatable::vel, max_pack, coll_prob
   double precision :: final_mass_array(NH), dds1, bin_frac1(NH), bin_frac2(NH)
   double precision, dimension(:), allocatable :: maxs, maxnc1, maxnc2, max_array1, maxim, timernd

   integer, dimension(:), allocatable :: check
   integer :: collindex, ncoll, pls
   integer :: f_counter, max_flag, bin11, bin22, f_flag
   double precision :: collselect(300), ptotal, prb_max,vvv
   double precision :: cmax1, cmax2, tcmax(NH)


if (rank1 .eq. 20) then
print*,dustmass,cell_dust_nden,cell_gas_den,cell_alpha,delta_time,coll_sch,coll_omega
end if

!  Set bin_frac as the fraction array in mass

   do i=1,NH
      bin_frac(i)=dustmass*in_arr(i)/(10d0**bin(i))
   end do

   bin_frac=bin_frac/sum(bin_frac)

!  dust_density=cell_dust_nden

   t = 0d0
   t2 = 0d0
   t3 = 0d0

   counter = 1
   tracker = 0

   f_flag = 0
   f_counter = 0
   max_flag = 0
   cmax1 = 0d0
   cmax2 = -9999999999d0

   bin11 = 0
   bin22 = 0


   do i=1,NH
      bin_wid_rs(i)=(10d0**bin_edge(i+1))-(10d0**bin_edge(i))
   end do

   max_array=0d0
    
   dust_num_den=cell_dust_nden

   dds1=cell_dust_nden


!  Initiate Kullback-Liebler divergence arrays: kbld1 and kbld2
!
!  We shall check the divergence after every 1000 time steps. If divergence is less than a 
!  pre-assigned value NSTAT declared in parameter.h file.
!
!  Cutting short the simulation using Kullback-Liebler Divergence will introduce error in the final 
!  distribution: For better precision, lower the value of NSTAT in parameter.h file (currently set at 0.1)


   do 
   
      tracker=tracker+1

      counter=counter+1

!  Select elements of non-zero fractions: Calculate number of particles to be drawn 
!  from each bin.

      non_zero_array=pack(bin_frac, (bin_frac .gt. 0d0))

      num_array=INT(bin_frac*dble(NTOTAL)+0.5d0)

!  Assign weights to each bin depending on their fraction:
!  If number of particles drawn from any bin < 1; We draw 1 particle from that bin.
!  The corresponding weight is set bu the original fraction.
!  For any bin with number of particles > 1; weight is always 1

      weights=0d0

      do i=1,NH

         if ((num_array(i) .eq. 0) .and. (bin_frac(i) .gt. 0d0)) then

            num_array(i)=1
            weights(i)=bin_frac(i)*dble(NTOTAL)

         elseif (num_array(i) .ge. 1) then

            weights(i)=1d0

         else

            weights(i)=0d0

         end if

      end do

!  Calculate dust number density for each bin. n_i = f_i x n_{dust}

      bin_num_den=bin_frac*dust_num_den

!if (rank1 .eq.20) then
!do nh1=22,30
!print*, bin_num_den(nh1),bin_frac(nh1),weights(nh1),num_array(nh1)
!end do
!end if

!  m_in_bin : Mass in the cell; same as the dust mass density of the cell
!             Not to mistake with mass in each bin. Same as M_{total} in paper
!             
!   *** This variable name will be changed later to avoid ambiguity ***
!
!            M_{total} = \Sum n_i x m_i
 
      m_in_bin=dot_product(bin_num_den,10d0**bin)

!  N = Total number of particles to be drawn from distribution
!
!  N depends on NTOTAL of parameter.h; NTOTAL-NH =< N =< NTOTAL+NH

      N=sum(num_array)

      allocate(array1(N))
      allocate(array2(N))
      allocate(par_frac_1(N))
      allocate(par_frac_2(N))
      allocate(weight1(N))
      allocate(weight2(N))
      allocate(coll_time(N))
      allocate(v(N))
      allocate(randarray(N))
      allocate(prob(N))


!  Populate arrays with particle from the distribution
!  
!  subroutine select_with_max: retains the maximum of sticking product from previous step

      call select_with_max(num_array, bin_edge, bin_frac, max_array, weights, N, dust_num_den, &
                           bin(2)-bin(1), array1, par_frac_1, weight1)
      
      call selection(num_array, bin_edge, bin_frac, weights, N, dust_num_den, &
                           bin(2)-bin(1), array2, par_frac_2, weight2)

       
!  Calculate relative velocity of collision: Ormel & Cuzzi 2007
!  
!  Also calculate probability of collision: p_{ij} = n_j x \sigma_{ij} x v_{ij}

      ptotal=0d0

      do nh1=1,NH
         if (bin_num_den(nh1) .gt. 0d0) then
            do nh2=1,NH
               if (bin_num_den(nh2) .gt. 0d0) then
                  
                  r1=(10d0**bin(nh1)/(FOURTHIRDS*PI*RHO_DUST))**ONETHIRD
                  r2=(10d0**bin(nh2)/(FOURTHIRDS*PI*RHO_DUST))**ONETHIRD

         ! Calculate collision cross-section: \sigma_{ij} = \pi (r_1 + r_2)^2
 
                  crosssection=PI*((r1+r2)**2d0)

                  call velocity(MIN(bin(nh1),bin(nh2)), MAX(bin(nh1),bin(nh2)), coll_rad, cell_z, &
                            cell_gas_den, cell_alpha, cell_temp, RHO_DUST, cell_cs, coll_sch, &
                            coll_omega, vvv)
!if (rank1 .eq. 20) then
!print*,'velocity',vvv,r1,r2!vvv,10d0**bin(nh1),10d0**bin(nh2)
!end if

                  ptotal=ptotal+(num_array(nh1)*num_array(nh2)*crosssection* &
                           vvv*dust_num_den/(2d0*dble(N)))
               end if
            end do
         end if
      end do

      do j=1,N

         ! Calculating particle sizes:

         r1=(10d0**array1(j)/(FOURTHIRDS*PI*RHO_DUST))**ONETHIRD
         r2=(10d0**array2(j)/(FOURTHIRDS*PI*RHO_DUST))**ONETHIRD

         ! Calculate collision cross-section: \sigma_{ij} = \pi (r_1 + r_2)^2
 
         crosssection=PI*((r1+r2)**2d0)

         call velocity(MIN(array1(j),array2(j)), MAX(array1(j),array2(j)), coll_rad, cell_z, &
                            cell_gas_den, cell_alpha, cell_temp, RHO_DUST, cell_cs, coll_sch, &
                            coll_omega, v(j))

         coll_time(j)=1d0/(par_frac_2(j)*crosssection*v(j))

         prob(j)=par_frac_2(j)*crosssection*v(j)

!        ptotal=ptotal+prob(j) !((1.26d-10)*crosssection*v(j)/((10d0**array2(j))*dble(N)))

      end do

!  Match O(n) total probability with the KMC O(n^2) total probability through scaling by N.

!     ptotal=ptotal*(dble(N))

      prb_max=maxval(prob) 

      deallocate(par_frac_1, par_frac_2)


      call randomnumber(randarray,N)

      randarray=randarray*prb_max

!  Collect pairs with successful collision and their weights, and relative velocity
!
!  Masses are followed in real space from this point
            
!  c1 & c2 : Successful collision pairs
!  wc1 & wc2 : Respective weights
!  vc : Relative velocity of collisions 

      c1=pack(10d0**array1, (randarray .lt. prob))
      c2=pack(10d0**array2, (randarray .lt. prob))

      wc1=pack(weight1, (randarray .lt. prob))
      wc2=pack(weight2, (randarray .lt. prob))

      vc=pack(v, (randarray .lt. prob))

!  nc1 & nc2 : Pairs that will not collide
!  wnc1 & wnc2 : Respective weights

      nc1=pack(10d0**array1, (randarray .gt. prob))
      nc2=pack(10d0**array2, (randarray .gt.prob))

      wnc1=pack(weight1, (randarray .gt. prob))
      wnc2=pack(weight2, (randarray .gt. prob))


      if (size(c1) .gt. NC) then

         coll_flag=1

         ncoll=NC

         allocate(check(size(c1)))

         check=0

         call randomnumber(collselect,300)

         do i=1,NC

            collindex=INT(collselect(i)*dble(size(c1)))+1

            check(collindex)=1

         end do

         s1=pack(c1, ((check .eq. 1) .and. (vc .lt. VFRAG)))
         s2=pack(c2, ((check .eq. 1) .and. (vc .lt. VFRAG)))
         ws1=pack(wc1, ((check .eq. 1) .and. (vc .lt. VFRAG)))
         ws2=pack(wc2, ((check .eq. 1) .and. (vc .lt. VFRAG)))

         
         f1=pack(c1, ((check .eq. 1) .and. (vc .gt. VFRAG)))
         f2=pack(c2, ((check .eq. 1) .and. (vc .gt. VFRAG)))
         wf1=pack(wc1, ((check .eq. 1) .and. (vc .gt. VFRAG)))
         wf2=pack(wc2, ((check .eq. 1) .and. (vc .gt. VFRAG)))

         nss1=pack(c1, (check .eq. 0))
         nss2=pack(c2, (check .eq. 0))
         wnss1=pack(wc1, (check .eq. 0))
         wnss2=pack(wc2, (check .eq. 0))

         deallocate(check)

      else

         coll_flag=0

         ncoll=size(c1)

         s1=pack(c1, (vc .lt. VFRAG))
         s2=pack(c2, (vc .lt. VFRAG))
         ws1=pack(wc1, (vc .lt. VFRAG))
         ws2=pack(wc2, (vc .lt. VFRAG))

         f1=pack(c1, (vc .gt. VFRAG))
         f2=pack(c2, (vc .gt. VFRAG))
         wf1=pack(wc1, (vc .gt. VFRAG))
         wf2=pack(wc2, (vc .gt. VFRAG))

      end if

     
      ! SET TIME EVOLUTION

         allocate(timernd(ncoll))

         call randomnumber(timernd, ncoll)

         timernd=log(timernd)

         t=t+((-1d0/ptotal)*sum(timernd))

         deallocate(timernd)

      
         bin11=INT((maxval([array1, array2])-bin_edge(1))/(bin_edge(2)-bin_edge(1)))+1

!if (rank1 .eq. 20) then         
! print*, maxval([array1,array2]),'  ',INT((maxval([array1, array2])-bin_edge(1))/(bin_edge(2)- &
!             bin_edge(1)))+1,'  ',t/SECYR, ptotal
!end if
!print*,MAXVAL([array1,array2]),' ',INT((maxval([array1,array2])-bin_edge(1))/(bin_edge(2)- &
!               bin_edge(1)))+1,'  ',t/SECYR,'  ',ncoll,'  ',tracker,'  ',size(ns1)

      w_from_f=0d0

      if (size(f1) .gt. 0) then

!        print*, 'Fragmentation ', size(f1)

         f_flag = 1         

         do i=1,size(f1)
            call new_frag(f1(i), wf1(i), f2(i), wf2(i), bin, bin_edge, &
                                    bin(2)-bin(1), temp_frag_mass)
            w_from_f=w_from_f+temp_frag_mass
         end do

      end if

      if (f_flag .eq. 1) then

         f_counter = f_counter+1

      end if


      ! COLLECT ALL WEIGHTS AND MARCH TOWARDS FINAL HISTOGRAM

      max_array=0d0

      if (coll_flag .eq. 0) then

         allocate(max_array1(size(s1)+size(nc1)+size(nc2)))
         max_array1=[s1+s2, nc1, nc2]

         do i=1,NH

            w_final(i)=0d0

            maxim=pack(max_array1, ((log10(max_array1) .gt. bin_edge(i)) .and. (log10(max_array1) .lt. &
                           bin_edge(i+1))))
            if (size(maxim) .gt. 0) then
               max_array(i)=maxval(maxim)
            end if

            hs=pack(s1*ws1+s2*ws2, ((log10(s1+s2) .gt. bin_edge(i)) .and. (log10(s1+s2) &
                      .lt. bin_edge(i+1))))
            if (size(hs) .gt. 0) then
               w_final(i)=w_final(i)+sum(hs)
            end if

            hns1=pack(nc1*wnc1, ((log10(nc1) .gt. bin_edge(i)) .and. (log10(nc1) .lt. bin_edge(i+1))))

            if (size(hns1) .gt. 0) then
               w_final(i)=w_final(i)+sum(hns1)
            end if

            hns2=pack(nc2*wnc2, ((log10(nc2) .gt. bin_edge(i)) .and. (log10(nc2) .lt. bin_edge(i+1))))

            if (size(hns2) .gt. 0) then
               w_final(i)=w_final(i)+sum(hns2)
            end if

            w_final(i)=w_final(i)+w_from_f(i)

            bin_num(i)=w_final(i)/(10d0**bin(i))

         end do

      else
         
         allocate(max_array1(size(s1)+2*size(nss1)+size(nc1)+size(nc2)))
         max_array1=[s1+s2, nss1, nss2, nc1, nc2]

         do i=1,NH
            
            w_final(i)=0d0
         
            maxim=pack(max_array1, ((log10(max_array1) .gt. bin_edge(i)) .and. (log10(max_array1) .lt. &
                           bin_edge(i+1))))
            if (size(maxim) .gt. 0) then
               max_array(i)=maxval(maxim)
            end if

            hs=pack(s1*ws1+s2*ws2, ((log10(s1+s2) .gt. bin_edge(i)) .and. (log10(s1+s2) &
                      .lt. bin_edge(i+1))))
            if (size(hs) .gt. 0) then
               w_final(i)=w_final(i)+sum(hs)
            end if

            hns1=pack(nss1*wnss1, ((log10(nss1) .gt. bin_edge(i)) .and. (log10(nss1) .lt. &
                     bin_edge(i+1))))

            if (size(hns1) .gt. 0) then
               w_final(i)=w_final(i)+sum(hns1)
            end if

            hns1=pack(nss2*wnss2, ((log10(nss2) .gt. bin_edge(i)) .and. (log10(nss2) .lt. &
                     bin_edge(i+1))))

            if (size(hns1) .gt. 0) then
               w_final(i)=w_final(i)+sum(hns1)
            end if

            hns1=pack(nc1*wnc1, ((log10(nc1) .gt. bin_edge(i)) .and. (log10(nc1) .lt. bin_edge(i+1))))

            if (size(hns1) .gt. 0) then
               w_final(i)=w_final(i)+sum(hns1)
            end if

            hns2=pack(nc2*wnc2, ((log10(nc2) .gt. bin_edge(i)) .and. (log10(nc2) .lt. bin_edge(i+1))))

            if (size(hns2) .gt. 0) then
               w_final(i)=w_final(i)+sum(hns2)
            end if

            w_final(i)=w_final(i)+w_from_f(i)

            bin_num(i)=w_final(i)/(10d0**bin(i))

         end do

      end if

      deallocate(max_array1)

      final_mass_array = (dustmass/sum(w_final))*w_final

      bin_frac1=0d0

      do i=1,NH

         bin_frac1(i) = final_mass_array(i)/(10d0**bin(i))

      end do

      bin_frac=(1d0/sum(bin_frac1))*bin_frac1


      dust_num_den=(dustmass/dot_product(10d0**bin, bin_frac))

      cmax1=final_mass_array(1)

      if (cmax1 .gt. cmax2) then

         max_flag = 1

         tcmax = final_mass_array

         cmax2 = cmax1
   
      end if

      if (bin22 .ne. bin11) then

         max_flag = 0

         f_counter = 0

         cmax2=-9999999999d0

      end if

      bin22 = bin11


      deallocate(array1, array2, v, prob)

      deallocate(weight1, weight2, coll_time, randarray)

!  Check condition for 1000 steps and call Kullback-Liebler divergence:

      if (f_counter .gt. nstep) then

         out_arr = tcmax

         dust_den_out = dust_num_den
      
         exit

      end if

   end do


end subroutine silicate_collision_SF






















































