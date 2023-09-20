#include "constants.h"
program name
    use vbtmodd
    implicit none
    real:: sigma,ohm,c_s,H,rhog,v_g,lamda,v_T,v_m,Re,t_L,t_ita
    real, dimension(500):: x,y
    real, dimension(500,500):: xm, ym, vbt, t1, t2, St1, St2, ep, ep1, vt_1, vt_2, vt
    integer:: i,j

    sigma=sigma0*(R/AU)**(-3/2)
    ohm=((G*MSUN)/R**3)**0.5
    c_s=(K_b*T/mump)**0.5
    H=c_s/ohm
    rhog=sigma/(H*(2*pi)**0.5)
    v_g=c_s*(alpha*8/pi)**0.5
    lamda=mump/(rhog*X_m)
    v_T=alpha*c_s*H
    v_m=0.5*lamda*c_s*(8/pi)**0.5
    Re=v_T/v_m
    t_L=1/ohm
    t_ita=t_L/(Re)**0.5
    call logspace(-4.0,4.0,500,x)
    call logspace(-4.0,4.0,500,y)
    call meshgrid(x,y,xm,ym)

    t1=(rhom*xm)/(rhog*c_g)
    t2=(rhom*ym)/(rhog*c_g)
    St1=t1/t_L
    St2=t2/t_L
    ep=St2/St1
    ep1=St1/St2

    do i = 1, 500
        do j = 1, 500
            if ( t1(i,j).LT.t_ita.AND.t1(i,j).GT.t2(i,j) ) then
                vt_1(i,j)=((v_g**2)*((St1(i,j)-St2(i,j))/(St1(i,j)+St2(i,j))))**0.5
                vt_2(i,j)=(((St1(i,j)**2)/(St1(i,j)+Re**(-0.5)))-(St2(i,j)**2)/(St2(i,j)+Re**(-0.5)))**0.5
                vt(i,j)=vt_1(i,j)*vt_2(i,j)
            else if (t2(i,j).LT.t_ita.AND.t2(i,j).GT.t1(i,j)) then
                vt_1(i,j)=((v_g**2)*((St2(i,j)-St1(i,j))/(St1(i,j)+St2(i,j))))**0.5
                vt_2(i,j)=(((St2(i,j)**2)/(St2(i,j)+Re**(-0.5)))-(St1(i,j)**2)/(St1(i,j)+Re**(-0.5)))**0.5
                vt(i,j)=vt_1(i,j)*vt_2(i,j)     
            else if (t_ita.LE.t1(i,j).AND.t1(i,j).LE.t_L.AND.t1(i,j).GT.t2(i,j)) then
                vt(i,j)=((v_g**2)*(St1(i,j))*(2*ya-1-ep(i,j)+(2/(1+ep(i,j)))* & 
                (1/(1+ya)+(ep(i,j)**3)/(ya+ep(i,j)))))**0.5
            else if (t_ita.LE.t2(i,j).AND.t2(i,j).LE.t_L.AND.t2(i,j).GT.t1(i,j)) then
                vt(i,j)=((v_g**2)*(St2(i,j))*(2*ya-1-ep1(i,j)+(2/(1+ep1(i,j)))* &
                (1/(1+ya)+(ep1(i,j)**3)/(ya+ep1(i,j)))))**0.5
            else if (t1(i,j).GT.t_L.AND.t1(i,j).GT.t2(i,j)) then
                vt(i,j)=((v_g**2)*(1/(1+St1(i,j))+1/(1+St2(i,j))))**0.5
            else if (t2(i,j).GT.t_L.AND.t2(i,j).GT.t1(i,j)) then 
                vt(i,j)=((v_g**2)*(1/(1+St2(i,j))+1/(1+St1(i,j))))**0.5           
            end if
        end do
    end do
    vbt=(vt**2+(8*K_b*T)*(rhom*(4/3)*pi*xm**3+rhom*(4/3)*pi*ym**3)/ &
    (pi*(rhom*(4/3)*pi*xm**3)*(rhom*(4/3)*pi*ym**3)))**0.5

    open(1, file='vbt.dat', status='unknown')
    do i = 1, 500
        write(1,*)vbt(:,i)
        
    end do
    close(1)
    open(1, file='xm.dat', status='unknown')
    do i = 1, 500
        write(1,*)xm(:,i)
        
    end do
    close(1)
    open(1, file='ym.dat', status='unknown')
    do i = 1, 500
        write(1,*)ym(:,i)
        
    end do
    close(1)

end program name