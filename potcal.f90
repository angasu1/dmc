      Module pot_calculation
      use types
      use global_variables   
     !use potentials
      Implicit none
      Contains

      subroutine potcalc(x,y,z,V)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        !Subroutine that calculates the potential for different molecules
        !it gets the cartesian coordinates for the cm-He distance in A. and
        !it returns the potential in cm-1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        real(rk),intent(in)::x,y,z
        real(rk),intent(out)::V
        real(rk)::xx(3)
 
          
       Select case (trim(molname))
         case('pyridine') 
         xx=(/x,y,z/)
         !v=PyridineHePot(xx)
         case('benzene')
         !call potential_benzene(x,y,z,v)
         case('ammonia') 
         !v=potw(x*ar2bo-xn(1),y*ar2bo-xn(2),z*ar2bo-xn(3),umb_ang)*har2cm
         !v=potw(x,y,z,umb_ang)  !OJO
         case('hcl')
         case('hf')
         case('hbr')
         case('hcn')
       end select



      end subroutine potcalc

              end Module pot_calculation

