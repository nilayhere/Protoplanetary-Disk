module vbtmodd
    implicit none
    
contains
    subroutine meshgrid(x,y,xm,ym)
    implicit none
    integer::i
    real, dimension(500):: x,y
    real, dimension(500,500):: xm,ym
    do i = 1, 500
        xm(:,i)=x
        ym(i,:)=y
    end do
    
    end subroutine meshgrid

    subroutine logspace(a,  b, n, x)
    implicit none
    real:: a,b, step
    real, dimension(n) :: x 
    integer :: i,n
    step=(abs(b-a))/(n-1)
    x(1)=10**(a)
    do i = 2, n
        x(i)=10**(a+(i-1)*step)
    end do
    
    end subroutine logspace
end module vbtmodd