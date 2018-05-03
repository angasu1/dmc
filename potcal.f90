    Module pot_calculation
      use types
      use global_variables   
      use pothx
      use math_subroutines
      use hcldimpot
      Implicit none

      Contains

      subroutine potcalc(y1,y2,y3,x,z,V)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        !Subroutine that calculates the potential for different molecules
        !it gets the cartesian coordinates for the centerofmass-He distance in A. and
        !it returns the potential in cm-1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        real(rk),intent(in)::x(3),z(3),y1(3),y2(3),y3(3)
        real(rk),intent(out)::V
        real(rk)::rvec(3),r,cthet
 
          
       Select case (trim(molname))
         case('pyridine') 
        !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hcl','hf','hbr','hcn')
         rvec=z-x
         r=dsqrt(dot_product(rvec,rvec))
         cthet=dot_product(rvec,y3)/r
         call potential(V,r,cthet)
       end select



      end subroutine potcalc

      subroutine potparam
       !Load the parameters for the corresponding potential

       Select case (trim(molname))
         case('pyridine') 
         !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hf')
         call pars_assignment(1)
         case('hcl')
         call pars_assignment(2)
         case('hbr')
         call pars_assignment(3)
         case('hcn')
       end select


      end subroutine potparam

      subroutine potdimparam
       !Load the parameters for the corresponding potential

       Select case (trim(molname))
         case('pyridine') 
         !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hf')
         case('hcl')
         call potkass
         case('hbr')
         case('hcn')
       end select


      end subroutine potdimparam

     subroutine potdimer(j,ri,rj,rvec,rad,potdim,xcm1,xcm2)
      integer(ik),intent(in)::j
      real(rk)::ri(3),rj(3),rvec(3),xcm1(3),xcm2(3)
      real(rk):: rad
      real(rk)::potdim,ctheta1,ctheta2,phi,cphi,temp1(3),temp2(3)
      real(rk)::cart_coor(3,4),fac1,fac2
      integer(ik)::i


     !rvec=(/3.746_rk,0.0_rk,0.0_rk/);rad=3.746_rk
     !ri=(/cos(9.0_rk*pi/180.0_rk),0.0_rk,sin(9.0_rk*pi/180.0_rk)/)
     !rj=(/cos(89.8_rk*pi/180.0_rk),0.0_rk,sin(89.8_rk*pi/180.0_rk)/)
     !xcm1=(/0.0_rk,0.0_rk,0.0_rk/)
     !xcm2=(/3.746_rk,0.0_rk,0.0_rk/)

      ctheta1=dot_product(rvec,ri)/rad
      ctheta2=dot_product(rvec,rj)/rad
      temp1=cross_unitary(ri,rvec)
      temp2=cross_unitary(rj,rvec)
      cphi=dot_product(-temp1,temp2)
      phi=dacos(cphi)
      if (abs(norm(ri)-1.0_rk).gt.1.0e-10_rk) stop 'norma1!!'
      if (abs(norm(rj)-1.0_rk).gt.1.0e-10_rk) stop 'norma2!!'
      fac1=-(att(2)%ma*requil/mtot)
      fac2=fac1+requil
     !write(*,*) 'fac',fac1,fac2
       cart_coor(:,1)=fac2*ri+xcm1
       cart_coor(:,2)=fac1*ri+xcm1
       cart_coor(:,3)=fac2*rj+xcm2
       cart_coor(:,4)=fac1*rj+xcm2
       if (j.eq.1) then
       call write_xyz(500,2,att%nom,cart_coor)
       endif 


       Select case (trim(molname))
         case('pyridine') 
         !xxvar=(/xvar,yvar,zvar/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hf')
         case('hcl')
           !ctheta1=dcos(9.0_rk*pi/180.0_rk);ctheta2=dcos(89.0_rk*pi/180.0_rk);phi=180.0_rk     
         call pothcl2(rad,ctheta1,ctheta2,phi,potdim)
           !write(*,*) 'rad',rad,ctheta1,ctheta2,phi,potdim
           !stop
         case('hbr')
         case('hcn')
       end select


     !write(*,*) 'pot',rad,ctheta1,ctheta2,phi,potdim
     !read*


      
           

            end subroutine potdimer

    end Module pot_calculation

