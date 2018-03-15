      program potencial
  integer, parameter :: rk = selected_real_kind(15,307)!real kind
  integer, parameter :: ik = selected_int_kind(8)!integer kind
  integer, parameter :: ark = selected_real_kind(15,307)!high-precision real kind
double precision r2,x,V,pi
Integer j,k
save j,k 
OPEN(UNIT=4,FILE='Potencial.txt')
pi=dacos(-1.0d0)
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
    
    
      subroutine potential(V,r2,x)
      implicit none
      integer l
      double precision  r2, x
      double precision fa, fb, Vsh, Vas, bohr_A,cm_H,c(6:7,1:3),R,S,T,F,suma,V,r1
!     
      Bohr_A=0.52917725d0
      cm_H=4.55634d-6
!     calculamos f_6 e f_7      
      fa=1-dexp(S(x)*(r2*Bohr_A))*suma(r2,x)
!     print*,'fa',fa
      fb=1-dexp(S(x)*(r2*Bohr_A))*(suma(r2,x)+((dabs(S(x)*(r2*Bohr_A)))**7)/5040) 
!     print*,'fb=',fb
      Vsh=T(r2,x)*dexp(R(x)+S(x)*(r2*Bohr_A))
!     print*,'Vsh=',Vsh
      Vas=(fa*c(6,1))/((r2*Bohr_A)**6)+(((fa*c(6,3))/((r2*Bohr_A)**6)))*F(2,x)
!     print*,'Vas=',Vas
      V=(Vsh+Vas)!*cm_H
!     hehfpot=(Vsh+Vas)*cm_H
      end
!     
!     P(n,x)
!     Polinomio de Legendre de orden n en x
!     
      double precision function F(l,x)
      implicit none
      integer l
      double precision x
      end
!     
!     Aqui calculamos e gardamos B(x), D(x) y G(x)
      double precision function R(x)
      implicit none
      integer y
      double precision  x,f
      R=0.d0
      do y = 1,6 
         R=R+b(y)*(F(y-1,x))
!     print*,'R(',x,')=',R
      enddo
      end
      
      double precision function S(x)
      implicit none
      double precision x,f
      integer l 
      S=0d0
      do l = 1,6
         S=S+D(l)*(F(l-1,x))
      enddo
      end

      double precision function T(r2,x)
      implicit none
      double precision x,r2,bohr_A,f
      integer l 
      T=0.0d0
      Bohr_A= 0.52917725d0
      do l = 1,6
         T=T+(g(1,l)+g(2,l)*(r2*Bohr_A)+g(3,l)*((r2*Bohr_A)**2)+g(4,l)&
     &        *((r2*Bohr_A)**3))*F(l-1,x)
!     print*,'T(',r2,',',x,')=',T
      enddo
      end
!     
!     Aqui calculamos suma
!     
      double precision function suma(r2,x)
      implicit none
      double precision r2,bohr_A,x,S
      integer k
      Bohr_A= 0.52917725d0
      suma=0d0
      do k = 1,7
         suma=suma+((dabs(S(x)*(r2*Bohr_A)))**(k-1))/fact(k-1)
!     print *, 'suma(',r2,',',x,')= ', suma
      enddo
      end



    
