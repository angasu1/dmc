Module pot_hx
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol        
    use types
    implicit none
   integer(ik), parameter::l=6 
   real(rk)::r, theta 
   real(rk), dimension(1:l):: d, b 
   real(rk), dimension (0:3)::c 
   real(rk), dimension (4,l)::g 
   real(rk), dimension (0:l-1)::F 
   Integer(ik) Factorial(0:7) 


    Contains

 subroutine pars_assignment(ind)

   Select case (ind)
     case (2)
  b=(/13.32010238_rk,0.34286792_rk,0.29816363_rk,0.29425822_rk,-0.00934684_rk,-0.03540720_rk/)
  d=(/-2.10089330_rk,0.01023381_rk,-0.12817719_rk,-0.01081582_rk,0.02043170_rk,0.00698361_rk/)
  g(1,:)=(/4.22760910_rk,-0.22438673_rk,1.58886495_rk,0.31882907_rk,0.39685627_rk,0.17984821_rk/)
  g(2,:)=(/-2.92257112_rk,0.25162254_rk,-1.21660466_rk,-0.11149131_rk,-0.26753161_rk,-0.08211105_rk/)
  g(3,:)=(/0.68212733_rk,-0.08245550_rk,0.31177458_rk,0.00000000_rk,0.05921384_rk,0.00955425_rk/)
  g(4,:)=(/-0.05598598_rk,0.00891702_rk,-0.02776035_rk,0.00250968_rk,-0.00434989_rk,0.000000000_rk/)
  c(6,3)=-0.00313099d0
  c(6,1)=-0.00566357d0
  fact=(/1_rk,1_rk,2_rk,6_rk ,24_rk,120_rk,720_rk/)
  F=(/1d0 ,x ,0.5d0*(3d0*(x**2)-1d0) ,0.5d0*(5d0*(x**3)-3d0*x) ,&
   &(1d0/8d0)*(35d0*(x**4)-30d0*(x**2)+3d0) ,(1d0/8d0)*(63d0*&
   &(x**5)-70d0*(x**3)+15d0*x)/)

   end select




 end subroutine pars_assignment







          end module pot_hx          
