module potinter                         
!Module that calculates the potential He-HX with X=F,Cl,Br,CN (depending on the 
!indicator 1,2,3,4 in the subroutine pothemol 
!interpolación V_lam*P_lam =V(R,theta)       
   use types                            
   use global_variables                 
   use input_values                     
   implicit none                        
                                        
   Contains                             
!Subroutine interpolacion(r0,th0,vlam, plam, nLinea, vsum)
   subroutine interpol(Vsum,r,cthet)    
       real(rk),intent(out)            :: Vsum
       real(rk),intent(in)             :: r,cthet
       integer(ik)                     :: i, j
       real(rk)                        :: th1, th2, r1, r2 !, th0, r0, vsum
       real(rk), allocatable           :: pl(:), vl(:)
       real(rk), allocatable           :: vlam(:,:), plam(:,:)
       nLinea = 1001                    
       allocate(pl(nColumna-1))         
       allocate(vl(nColumna-1))         
                                        
       !interpolación, encontrar los vlam, plam y sacar el potencial, para nuevos r0 y th0
       do i = 1, nLinea                 
           if((plam(i,1) <= cthet).and.(cthet <= plam(i+1,1))) then
              th1 = plam(i,1)           
              th2 = plam(i+1,1)         
              do j = 2, nColumna        
                  pl(j-1) = plam(i,j) + ((cthet-th1)/(th2-th1))*(plam(i+1,j) - plam(i,j))
              end do !j = 2, nColumna          
           end if!(pl(i,1).lt.th0.gt.pl(i+1,1)) then
       end do ! i = 1, nlinea                  
       do i = 1, nLinea                 
           if((vlam(i,1) <= r).and.(r <= vlam(i+1,1))) then
              r1 = vlam(i,1)            
              r2 = vlam(i+1,1)          
              do j = 2, nColumna        
                 vl(j-1) = vlam(i,j) + ((r-r1)/(r2-r1))*(vlam(i+1,j) - vlam(i,j))
              end do !j = 2, nColumna          
              Vsum = sum(vl*pl)         
              ! write(*,*) vl(i), pl(i), vsum   
              ! read *                          
             end if!(pl(i,1).lt.th0.gt.pl(i+1,1)) then
       end do ! i = 1, nlinea                  
                                        
                                        
1000 format(*(f28.16,2x))               
                                        
   end subroutine interpol              
                                        
end module potinter     
