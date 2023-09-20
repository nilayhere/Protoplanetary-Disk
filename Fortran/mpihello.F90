program hellompi
    use mpi
    implicit none
    !include 'mpif.h'

    integer :: ierr,rrank,nprocs

    call MPI_INIT(ierr)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rrank, ierr)

    print *, "hello world I am processor ", rrank, "of", nprocs, "processors"

    call MPI_FINALIZE(ierr)

    
end program hellompi