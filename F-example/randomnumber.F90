!-----------------------------------------------------------------
!
!                         RANDOM NUMBER GENERATOR
!
!  subroutine  randomnumber(r,n)
!              r: array of random numbers
!               n: dimension of array r
!-----------------------------------------------------------------


   subroutine randomnumber(r,n)

        integer, intent(in) :: n
        double precision :: r(n)
        integer :: i,j,pp
        double precision :: rand
        call init_random_seed()
        call random_number(r)
   
   end subroutine randomnumber


         subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            !integer, parameter :: k=selected_int_kind(16)
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
            integer :: getpid
          
            call random_seed(size = n)
            allocate(seed(n))
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed




 !------------------------------------------------------------------------------------------
