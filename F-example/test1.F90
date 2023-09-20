program name
    implicit none
    real, allocatable, dimension(:):: x
    allocate(x(5))
    x(1)=6
    print*, x
end program name