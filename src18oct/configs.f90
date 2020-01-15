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

        function radius(a,b)
        real(rk)::a(3),b(3),c(3),radius
        c=b-a
        radius=dsqrt(dot_product(c,c))
        return
        end function radius

        function deg(rad)
         real(rk)::deg,rad

           deg=rad*180_rk/pi

          return
         end function deg       

!        function norm(a)
!        real(rk)::norm,a(3)
!        radius=dsqrt(dot_product(a,a))

!       return
!       end function norm

             subroutine rotz(rij,th)
               real(rk)::rij(3),th,Amat(3,3),rijt(3)
               integer(ik)::i,j

               amat(1,:)=(/dcos(th),dsin(th),0.0_rk/)
               amat(2,:)=(/-dsin(th),dcos(th),0.0_rk/)
               amat(3,:)=(/0.0_rk,0.0_rk,1.0_rk/)
               rijt=0.0_rk

               do i=1,3
                       do j=1,3
                         rijt(i)=rijt(i)+amat(i,j)*rij(j)
                       enddo
               enddo

                  rij=rijt      
                               

              return
              end subroutine rotz      

             subroutine roty(rij,th)
               real(rk)::rij(3),th,Amat(3,3),rijt(3)
               integer(ik)::i,j

               amat(1,:)=(/dcos(th),0.0_rk,dsin(th)/)
               amat(2,:)=(/0.0_rk,1.0_rk,0.0_rk/)
               amat(3,:)=(/-dsin(th),0.0_rk,dcos(th)/)
               rijt=0.0_rk

               do i=1,3
                       do j=1,3
                         rijt(i)=rijt(i)+amat(i,j)*rij(j)
                       enddo
               enddo

                  rij=rijt      
                               

              return
              end subroutine roty      


end module configur

  program confs
   use configur
   implicit none
   integer(ik)::i,j,k,idum
   real(rk)::y(3,2),xx(3,4),rij(3),theta(2),phi(2),signo,xcm(3,2)
   real(rk)::m1,m2,mtot,requil,fac1,fac2,boxl,throt,phrot
   real(rk)::ctheta1,ctheta2,temp1(3),temp2(3),cphi,phi1,phi2,phi3,phi4,error
   open(unit=100,file='cart.xyz')
   open(unit=200,file='angles.dat')
   open(unit=300,file='radius.dat')

   rij=(/3.8_rk,0.0_rk,0.0_rk/)
   m1=1.0_rk
   m2=35.453_rk
   mtot=m1+m2
   requil=1.5_rk
   idum=-12345
   fac1=-(m1*requil/mtot)
   fac2=fac1+requil
   boxl=10.0_rk

configs:do i=1,100
            !Calculation of rij vector
            do j=1,2
                    do k=1,3
                      xcm(k,j)=(0.5_rk-ran2(idum))*boxl      
                    enddo
            enddo      

             rij=xcm(:,2)-xcm(:,1)
             
              

     !Calculation of yi vecs
     do j=1,2      
      theta(j)=ran2(idum)*pi
      phi(j)=ran2(idum)*2.0_rk*pi
      y(:,j)=(/dsin(theta(j))*dcos(phi(j)),dsin(theta(j))*dsin(phi(j)),dcos(theta(j))/)  
     enddo

     !Rotation to the system
             phrot=datan2(-rij(2),rij(1))
             throt=dacos(rij(3)/norm(rij))
             throt=pi/2.0_rk-throt

            !write(*,*) 'rij ant',rij
             call rotz(rij,-phrot)
             call rotz(y(:,1),-phrot)
             call rotz(y(:,2),-phrot)
            !write(*,*) 'rij desp',rij
             call roty(rij,throt)
             call roty(y(:,1),throt)
             call roty(y(:,2),throt)
            !write(*,*) 'rij desp2',rij
            !stop



     !Calculation of cartesian coords
       xx(:,1)=fac2*y(:,1)-rij/2.0_rk
       xx(:,2)=fac1*y(:,1)-rij/2.0_rk
       xx(:,3)=fac2*y(:,2)+rij/2.0_rk
       xx(:,4)=fac1*y(:,2)+rij/2.0_rk


     !Angles for potential calculation
      ctheta1=dot_product(y(:,1),rij)/norm(rij)
      ctheta2=dot_product(y(:,2),rij)/norm(rij)

      !phi, method 1
      temp1=cross_unitary(y(:,1),rij)
      temp2=cross_unitary(y(:,2),rij)
      cphi=dot_product(temp1,temp2)
      phi4=dacos(cphi)

      !Phi method 2
      phi1=datan2(xx(2,1),xx(3,1))
      phi2=datan2(xx(2,3),xx(3,3))
      phi3=abs(phi1-phi2)

      if(abs(deg(phi3)-deg(phi4)).lt.1.0e-10_rk) then
        error=0.0_rk
      else
        error=abs(deg(phi3)-360.0_rk+deg(phi4))
      endif


      !Writing in files
      call write_xyz(100,4,atnom,xx)
      write(200,1000) deg(acos(ctheta1)),deg(acos(ctheta2)),deg(phi4),deg(phi3),error
      write(300,1000) radius(xx(:,1),xx(:,2)),radius(xx(:,3),xx(:,4))&
                      &,radius((m1*xx(:,1)+m2*xx(:,2))/mtot,(m1*xx(:,3)+m2*xx(:,4))/mtot)


 1000 format(*(F18.10,3x))
    
      enddo configs


   end program confs        
