module potinter
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol 
!interpolación V_lam*P_lam =V(R,theta)       
   use types
   use global_variables 
   use input_values
   implicit none

   Contains
   subroutine interpol(Vsum,r,cthet)
       real(rk),intent(out)            :: Vsum
       real(rk),intent(in)             :: r,cthet
       integer(ik)                     :: i, j
       real(rk)                        :: cth1, cth2, r1, r2 !, th0, r0, vsum
       integer(ik)                     :: linear, lineath
      !real(rk), allocatable           :: pl(:), vl(:)
       real(rk)                        :: pl, vl,slopeth,sloper
       
       if (r.gt.rexit) then
               vsum=0.0_rk
               return
       endif



       nColumna = ntermspot + 2
      !allocate(pl(nColumna-1))                
      !allocate(vl(nColumna-1))                
       !interpolación, encontrar los vlam, plam y sacar el potencial, para nuevos r0 y th0
       linear = 1 + int((1000._rk*(r-2._rk))/38._rk)
       !lineath = 1 + int((1000._rk*(cthet))/pi) 
       lineath = 1 - nint(500._rk*(cthet-1._rk)) 
       

       cth1 = plam(lineath,1)                
       cth2 = plam(lineath+1,1)              
       r1 = vlam(linear,1)                  
       r2 = vlam(linear+1,1)                

       ! Asegura estar en la linea correcta 
       if(abs(cthet-cth1)>(1._rk/1000._rk).or.abs(r-r1)>(38._rk/1000._rk))then
           write(*,*) "no estas en la línea correcta"
           stop                                 
       end if         

       if (lineath.gt.1001) stop 'fuera de rango theta'
       if (linear.gt.1001) stop 'fuera de rango r'

       ! encuentra valores nuevos de vlam y plam
       vsum=0
       slopeth=(cthet-cth1)/(cth2-cth1)
       sloper=(r-r1)/(r2-r1)
       do j = 2, nColumna                   
           pl = plam(lineath,j) + slopeth*(plam(lineath+1,j) - plam(lineath,j))
           vl = vlam(linear,j)  + sloper*(vlam(linear+1,j) - vlam(linear,j))
           vsum=vsum+pl*vl
       end do !j = 2, nColumna              
                                        
        write(*,*) 'pot=',vsum,r,cthet
        read*
                                        
   10 format(*(f28.16,2x))               

   end subroutine interpol        

end module potinter
