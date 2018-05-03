module configur
implicit none

  integer, parameter :: rk = selected_real_kind(15,307)!real kind
  integer, parameter :: ik = selected_int_kind(8)!integer kind
  integer, parameter :: ark = selected_real_kind(15,307)!high-precision real kind

character(len=2)::atnom(4)=(/' H','Cl',' H','Cl'/)
integer(ik)::nat=4
real(rk)::pi=dacos(-1.0_rk)
Contains


       FUNCTION ran2(idum) 
       integer,intent(inout)::idum        
       integer(ik) IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
       Real (rk) ran2,AM,EPS,RNMX 
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1._rk/IM1,IMM1=IM1-1,& 
      &IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,& 
      &NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7_rk,RNMX=1._rk-EPS) 
       integer(ik) idum2,j,k,iv(NTAB),iy 
       SAVE iv,iy,idum2 
       DATA idum2/123456789/, iv/NTAB*0/, iy/0/ 
       if (idum.le.0) then 
         idum=max(-idum,1) 
         idum2=idum 
         do 11 j=NTAB+8,1,-1 
  
           k=idum/IQ1 
           idum=IA1*(idum-k*IQ1)-k*IR1 
           if (idum.lt.0) idum=idum+IM1 
           if (j.le.NTAB) iv(j)=idum 
 11      continue 
         iy=iv(1) 
       endif 
       k=idum/IQ1 
       idum=IA1*(idum-k*IQ1)-k*IR1 
       if (idum.lt.0) idum=idum+IM1 
       k=idum2/IQ2 
       idum2=IA2*(idum2-k*IQ2)-k*IR2 
       if (idum2.lt.0) idum2=idum2+IM2 
       j=1+iy/NDIV 
       iy=iv(j)-idum2 
       iv(j)=idum 
       if(iy.lt.1)iy=iy+IMM1 
       ran2=min(AM*iy,RNMX) 
       return 
       END FUNCTION ran2 

    function cross_unitary(a,b)
    real (rk) ::a(3),b(3),cross_unitary(3),anorm
    cross_unitary=0._rk

    cross_unitary(1)=a(2)*b(3)-a(3)*b(2)
    cross_unitary(2)=a(3)*b(1)-a(1)*b(3)
    cross_unitary(3)=a(1)*b(2)-a(2)*b(1)
    anorm=norm(cross_unitary)
    if (anorm.lt.1.e-20_rk) stop 'problems calculating cross_unitary'
    cross_unitary=cross_unitary/anorm
    return
    end function cross_unitary
     function norm(a)
             real(rk)::a(3),norm

             norm=dsqrt(dot_product(a,a))


             end function norm

             subroutine write_xyz(nfile,nat,atyp,xx)
          integer(ik)::i,j,nfile,nat      
          real(rk)::xx(3,2*nat)
          character(len=2)::atyp(nat)
 
          write(nfile,*) nat
          write(nfile,*) 
 
         do i=1,nat
             write(nfile,1000) atyp(i),(xx(j,i),j=1,3)
         enddo


 1000    format(A,3(F20.6))
         return
        end subroutine write_xyz


end module configur

  program confs
   use configur
   implicit none
   integer(ik)::i,j,idum
   real(rk)::y(3,2),xx(3,4),rij(3),theta(2),phi,signo
   real(rk)::m1,m2,mtot,requil,fac1,fac2
   real(rk)::ctheta1,ctheta2,temp1(3),temp2(3),cphi
   open(unit=100,file='cart.xyz')
   open(unit=200,file='angles.dat')

   rij=(/3.8_rk,0.0_rk,0.0_rk/)
   m2=35.453_rk
   m1=1.0_rk
   mtot=m1+m2
   requil=1.5_rk
   idum=-12345
   fac1=-(m1*requil/mtot)
   fac2=fac1+requil

    do i=1,100
     do j=1,2      
      theta(j)=ran2(idum)*pi
      signo=1.0_rk      
      if (ran2(idum).lt.0.5_rk) signo=-1.0_rk
      y(:,j)=(/dsin(signo*theta(j)),dcos(signo*theta(j)),0.0_rk/)  
     !write(*,*) 'fac',fac1,fac2
      if (j.eq.1) then
       xx(:,1)=fac2*y(:,j)
       xx(:,2)=fac1*y(:,j)
      else
       xx(:,3)=fac2*y(:,j)+rij
       xx(:,4)=fac1*y(:,j)+rij
      endif
     enddo
      ctheta1=dot_product(y(:,1),rij)/3.8_rk
      ctheta2=dot_product(y(:,2),rij)/3.8_rk
      temp1=cross_unitary(y(:,1),rij)
      temp2=cross_unitary(y(:,1),rij)
      cphi=dot_product(-temp1,temp2)
      phi=dacos(cphi)
             call write_xyz(100,4,atnom,xx)
             write(200,1000) dcos(theta(1))-ctheta1,dcos(theta(2))-ctheta2,phi*180_rk/pi
 1000 format(*(F18.10,3x))
    enddo


   end program confs        
