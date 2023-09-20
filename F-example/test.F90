program name
    implicit none
    type :: zdata
      
    double precision :: n
    double precision :: rho_g, temperature, alpha, c_s, zvalue, rho_d
    double precision :: hist, delz
    end type zdata

    type(zdata), dimension(:), allocatable , target :: z
    allocate(z(2))
    z(1)%n=5
    z(2)%n=7
    READ*, z(1)%hist
    print*, z
    print*, z(1)%hist

end program name