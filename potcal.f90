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

     subroutine potdimer(ri,rj,rvec,rad,potdim)
      real(rk),intent(in)::ri(3),rj(3),rvec(3)
      real(rk),intent(in):: rad
      real(rk)::potdim,ctheta1,ctheta2,phi,cphi,temp1(3),temp2(3)


      ctheta1=dot_product(rvec,ri)/rad
      ctheta2=dot_product(rvec,rj)/rad
      temp1=cross_unitary(rvec,ri)
      temp2=cross_unitary(-rvec,rj)
      cphi=dot_product(temp1,temp2)
      phi=dacos(cphi)

      call pothcl2(rad,ctheta1,ctheta2,phi,potdim)

     !write(*,*) 'pot',rad,ctheta1,ctheta2,phi,potdim
     !read*


      
           

            end subroutine potdimer

    end Module pot_calculation

