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
       real(rk)                        :: th1, th2, r1, r2 !, th0, r0, vsum
       integer(ik)                     :: linear, lineath
       real(rk), allocatable           :: pl(:), vl(:)
       
       nColumna = ntermspot + 2
       allocate(pl(nColumna-1))                
       allocate(vl(nColumna-1))                
                                 
       !interpolación, encontrar los vlam, plam y sacar el potencial, para nuevos r0 y th0
       linear = 1 + int((1000._rk*(r-2._rk))/33._rk)
       !lineath = 1 + int((1000._rk*(cthet))/pi) 
       lineath = 1 - nint(500._rk*(cthet-1._rk)) 
       

       th1 = plam(lineath,1)                
       th2 = plam(lineath+1,1)              
       r1 = vlam(linear,1)                  
       r2 = vlam(linear+1,1)                

       ! Asegura estar en la linea correcta 
       if(abs(cthet-th1)>(1._rk/1000._rk).or.abs(r-r1)>(33._rk/1000._rk))then
           write(*,*) "no estas en la línea correcta"
           stop                                 
       end if         

       ! encuentra valores nuevos de vlam y plam
       do j = 2, nColumna                   
           pl(j-1) = plam(lineath,j) + ((cthet-th1)/(th2-th1))*(plam(lineath+1,j) - plam(lineath,j))
       end do !j = 2, nColumna              
                                        
       do j = 2, nColumna                   
           vl(j-1) = vlam(linear,j) + ((r-r1)/(r2-r1))*(vlam(linear+1,j) - vlam(linear,j))
       end do !j = 2, nColumna 
       vsum = sum(vl*pl)                
                                        
                                        
   10 format(*(f28.16,2x))               

   end subroutine interpol        

end module potinter
