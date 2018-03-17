    Module pot_calculation
      use types
      use global_variables   
      use pothx
      Implicit none

      Contains

      subroutine potcalc(y1,y2,y3,x,z,V)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        !Subroutine that calculates the potential for different molecules
        !it gets the cartesian coordinates for the cm-He distance in A. and
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
         case('hcl')
         rvec=z-x
         r=dsqrt(dot_product(rvec,rvec))
         cthet=dot_product(rvec,y3)/r
         call potential(V,r,cthet)
         case('hf')
         rvec=z-x
         r=dsqrt(dot_product(rvec,rvec))
         cthet=dot_product(rvec,y3)/r
         call potential(V,r,cthet)
         case('hbr')
         rvec=z-x
         r=dsqrt(dot_product(rvec,rvec))
         cthet=dot_product(rvec,y3)/r
         call potential(V,r,cthet)
         case('hcn')
       end select



      end subroutine potcalc

      subroutine potparam

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

    end Module pot_calculation

