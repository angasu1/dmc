module potinter
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol        
   use types
   use input_values
   implicit none



    Contains


      subroutine interpol(V,r,cthet)
        real(rk),intent(out)::V
        real(rk),intent(in)::r,cthet




      end subroutine interpol        




   end module potinter
