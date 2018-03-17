Module pot_hx
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol        
    implicit none
  integer, parameter :: rk = selected_real_kind(15,307)!real kind
  integer, parameter :: ik = selected_int_kind(8)!integer kind
   integer(ik), parameter::l=5 
   real(rk)::r, theta 
   real(rk), dimension(0:l):: d, b 
   real(rk), dimension (2,2)::c 
   real(rk), dimension (0:3,0:l)::g 
   real(rk), dimension (0:l)::F 
   real(rk),parameter::bo2ar=0.5291772_rk
   Integer(ik) fact(0:7) 


    Contains

 subroutine pars_assignment(ind)
 integer(ik)::ind

   Select case (ind)
     case (2)
  b=(/13.32010238_rk,0.34286792_rk,0.29816363_rk,0.29425822_rk,-0.00934684_rk,-0.03540720_rk/)
  d=(/-2.10089330_rk,0.01023381_rk,-0.12817719_rk,-0.01081582_rk,0.02043170_rk,0.00698361_rk/)
  g(0,:)=(/4.22760910_rk,-0.22438673_rk,1.58886495_rk,0.31882907_rk,0.39685627_rk,0.17984821_rk/)
  g(1,:)=(/-2.92257112_rk,0.25162254_rk,-1.21660466_rk,-0.11149131_rk,-0.26753161_rk,-0.08211105_rk/)
  g(2,:)=(/0.68212733_rk,-0.08245550_rk,0.31177458_rk,0.00000000_rk,0.05921384_rk,0.00955425_rk/)
  g(3,:)=(/-0.05598598_rk,0.00891702_rk,-0.02776035_rk,0.00250968_rk,-0.00434989_rk,0.000000000_rk/)
  c(1,1)=-0.00566357d0;c(1,2)=-0.00313099d0;c(2,1)=0.0_rk;c(2,2)=0.0_rk
  fact=(/1_rk,1_rk,2_rk,6_rk ,24_rk,120_rk,720_rk,5040_rk/)

   end select
   end subroutine pars_assignment

      subroutine potential(V,r,x)
      real(rk),intent(in)::  r, x
      real(rk):: fn(2), Vsh, Vas, bohr_A,cm_H,V
      real(rk)::x2,x3,x4,x5,bb,dd,gg,xf,ra,r2,r6,r7
      integer(ik)::i
!     
      ra=r*bo2ar
      x2=x*x;x3=x2*x;x4=x2*x2;x5=x4*x
      F=(/1._rk ,x ,0.5_rk*(3._rk*(x2)-1_rk) ,0.5_rk*(5._rk*(x3)-3._rk*x) ,&
      &(1._rk/8._rk)*(35._rk*(x4)-30._rk*(x2)+3._rk) ,(1._rk/8._rk)*(63._rk*&
      &(x5)-70._rk*(x3)+15._rk*x)/)

      bb=bfun()
      dd=dfun()
      gg=gfun(ra)
      Vsh=gg*dexp(bb+dd*ra)
      xf=dd*ra
      fn=f67(xf)

      Vas=0.0_rk
      r2=ra*ra;r6=r2*r2*r2;r7=r6*ra
      Vas=fn(1)*(c(1,1)+c(1,2)*F(2))/R6       
      Vas=vas+fn(2)*(c(2,1)*F(1)+c(2,2)*F(3))/R7       
      V=(Vsh+Vas)!*cm_H
     !write(*,*) 'pots',fn
     !write(*,*) 'pots',Vas,Vsh,V
     !stop
!     hehfpot=(Vsh+Vas)*cm_H
      end

      function bfun()
        real(rk)::bfun
        integer(ik)::k

        bfun=0.0_rk
        do k=0,l
          bfun=bfun+b(k)*F(k)
        enddo

        return
       end function bfun

      function dfun()
        real(rk)::dfun
        integer(ik)::k

        dfun=0.0_rk
        do k=0,l
          dfun=dfun+d(k)*F(k)
        enddo

        return
       end function dfun

      function gfun(r)
        real(rk)::gfun,r,r2,r3
        integer(ik)::k

        gfun=0.0_rk
        r2=r*r;r3=r2*r
        do k=0,l
          gfun=gfun+(g(0,k)+g(1,k)*r+g(2,k)*r2+g(3,k)*r3)*f(k)
        enddo

        return
       end function gfun

       function f67(x)
        real(rk)::f67(2),x,abx,x2,x3,x4,x5,x6,x7,suma(2)
        integer(ik)::l

        abx=abs(x);x2=abx*abx;x3=abx*x2;x4=x2*x2;x5=x4*abx;x6=x4*x2;x7=x4*x3
        suma(1)=1.0_rk+x+x2/fact(2)+x3/fact(3)+x4/fact(4)+x5/fact(5)+x6/fact(6)
        suma(2)=suma(1)+x7/fact(7)
        suma=suma*dexp(x)      
        f67=1.0_rk+suma      


        return
       end function f67

       end module pot_hx          







program potencial
use pot_hx
double precision r2,x,V,pi
Integer j,k
save j,k 
OPEN(UNIT=4,FILE='Potencial.txt')
pi=dacos(-1.0d0)
call pars_assignment(2)
do j=0,100
   write(*,*) j
r2=6d0+j/20d0
do i=0,100
x=1.0d0-i/50.0d0
call potential(V,r2,x)
WRITE(4,50) dacos(x)*180_rk/pi,r2,V
      50     FORMAT(F18.10,3X,F18.10,3X,F18.10,3X)

End do
      !write(4,*)
end do

End
    



    
