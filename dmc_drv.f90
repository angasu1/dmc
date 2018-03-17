Module dmc_cal 
  use types
  use global_variables      
  use leg_integrator
  use input_values 
  use math_subroutines
  use phys_cons
  use pot_calculation 
  implicit none      
 Save
 integer,parameter::delta_N_max=200
 real(rk),allocatable:: y1(:,:),y2(:,:),y3(:,:)
 real(rk),allocatable:: zcor(:,:,:),xcor(:,:)
 real(rk),allocatable::Er(:),pot(:),eloc(:),erme_wn(:)
 real(rk),allocatable::signo(:),fonda(:,:),alphahe(:),betahe(:)
 real(rk)::ermedia,ermedia2
 real(rk),parameter::cc=20000d0,aa=0.5d0
 real(rk)::alpham,betam,gammam
 real(rk)::alphahe1,betahe1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!Programa que usa el metodo de DMC para calcular estados cuanticos! 
!y la funcion de onda correspondiente para el complejo He_NH3     !
!!!!!!!!!!!!!!,xx,mass!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

Contains

Subroutine dmc_drv 
  real(rk)    :: media,sd
  integer(ik) :: icorrida,N,npaso,j,icontr

   allocate(y1(3,nw+delta_N_max),y2(3,nw+delta_N_max),y3(3,nw+delta_N_max))
   allocate(zcor(3,nw+delta_N_max,nhe),xcor(3,nw+delta_N_max))
   allocate(Er(0:nstps),pot(nw+delta_N_max),eloc(nw+delta_N_max))
   allocate(signo(nw+delta_N_max),fonda(nw+delta_N_max,0:1))
   allocate(alphahe(nhe),betahe(nhe),erme_wn(nruns))


       idum=-1234
       if (.not.fixran) call seed_cal(idum)
       call rotcons_calculation
       call potparam

 main_loop:  do icorrida=1,nruns
     

   !initialization of values
         icontr=0               
         ermedia=0d0
         ermedia2=0d0 
         erme_wn(icorrida)=0d0
         npaso=1 
         er(0)=0d0 
         N=nw  
       
        call initial_conditions(icorrida) 

       !Ciclo sobre todos los pasos temporales 

steps:  do npaso=1,nstps
       !Se llama a la subrutina que efectúa la difusión
        call diffusion(npaso,N)
       !Se calcula el potencial de cada replica
        pot=0d0
        call potential_calculation(npaso,N,icontr) 
!       !Se llama a la subrutina que efectúa los procesos de vida o muerte
        call branch(npaso,N,icorrida) 
        enddo steps


!       !Se escriben los valores de la energía media final en cada corrida
        write (100,*) erme_wn(icorrida)



enddo main_loop

        call write_emedia
        Close(100)
        close(250)
        Close(300)
        close(400)   

      


End Subroutine dmc_drv

      Subroutine initial_conditions(icorrida)
      real (rk)   :: rr,u,phiinicial,rot,rmin,rmax
      real (rk)   :: phi,theta,psi,rrot(3),rhe(3,nw+delta_N_max,nhe)
      integer :: jhe,j,icorrida,icont
 
    icont=0
    !Primero. Se asignan las posiciones del centro de masas de la molecula 
        xcor=0D0 
        

  !Segundo. Se asignan los vectores de orientacion de la molecula 
  
       do j=1,nw 
       y1(:,j)=ei
       y2(:,j)=ej
       y3(:,j)=ek

      rot=1.0d0 !poner a cero si se quiere iniciar alineado con el lab
      phiinicial = ran2(idum)*2.d0*pi*rot
       !Se rota la molecula con respecto al ejex molecular un angulo phiinicial 
      call rotacion(y1(1,j),Y1(2,j),Y1(3,j),y3(1,j),y3(2,j),y3(3,j),phiinicial) 
      call rotacion(y1(1,j),Y1(2,j),Y1(3,j),y2(1,j),y2(2,j),y2(3,j),phiinicial) 
        
        phiinicial = ran2(idum)*2.d0*pi*rot
       !se rota la molecula con respecto al ejey molecular un angulo phiinicial 
      call rotacion(y2(1,j),Y2(2,j),Y2(3,j),y3(1,j),y3(2,j),y3(3,j),phiinicial) 
      call rotacion(y2(1,j),Y2(2,j),Y2(3,j),y1(1,j),y1(2,j),y1(3,j),phiinicial) 
        

      if (moltyp.gt.2) then !Rotate for a non linear molecule
       phiinicial = ran2(idum)*2.d0*pi*rot 
       !se rota la molecula con respecto al ejez molecular un angulo phiinicial 
       call rotacion(y3(1,j),Y3(2,j),Y3(3,j),y1(1,j),y1(2,j),y1(3,j),phiinicial) 
       call rotacion(y3(1,j),Y3(2,j),Y3(3,j),y2(1,j),y2(2,j),y2(3,j),phiinicial)
      endif

     call euler(y1(:,j),y2(:,j),y3(:,j),phi,theta,psi)
 !      
        
!       !Tercero. Se asignan los valores de las posiciones iniciales del helio aleatoriamente
!       !en una esfera de radio rr cerca del mínimo de potencial  
          do jhe=1,nhe
            rmin=6.0_rk
            rmax=15.0_rk
            if (simtyp.eq.2) then
            rmin=rfb;rmax=rfb
            endif       
            rr=rmin+ran2(idum)*(rmax-rmin)        
            u = ran2(idum)*2.d0-1.d0 
            phiinicial = ran2(idum)*2.d0*pi 
                             
            zcor(1,j,jhe) = rr*dsqrt(1-u**2)*dcos(phiinicial) 
            zcor(2,j,jhe) = rr*dsqrt(1-u**2)*dsin(phiinicial) 
            zcor(3,j,jhe) = rr*u  
      
           rhe(:,j,jhe)=(zcor(:,j,jhe)-xcor(:,j))
       
           call rotaciones(rhe(:,j,jhe),rrot,phi,theta,psi)
           
            if (icont.lt.10000) then 
             write (400,1000) rrot,rr,phi,theta,psi
             icont=icont+1
            endif
          enddo 
             enddo
            write(400,*) 'New run'

1000 format(7(F10.4))

        return
     end Subroutine initial_conditions
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!Esta subrutina hace dar un paso a todos los caminantes, rota la mo-!
!lécula,  manda a calcular el potencial y recalcula                 ! 
!sus posiciones                                                     ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     Subroutine diffusion(npaso,N) 
       real (rk),dimension(3)::Fq_x 
       real (rk),dimension(3,nhe)::Fq_z 
       real (rk)  :: phix,phiy,phiz,rvec_he(3)
       integer     :: npaso,j,i,jhe
       integer,intent(in)::N 
!        
       do j=1,N 
!        
        If (is) Then
        Call Fq(j,Fq_z,Fq_x)
        Else
        Fq_z=0d0
        Fq_x=0d0
        endif
        
        if (simtyp.eq.2) then
                do jhe=1,nhe 
                rvec_he=zcor(:,j,jhe)-xcor(:,j)        
       Phix=dsqrt(2d0*dtau*diffcHe)*gasdev(idum)!*0d0 
       Phiy=dsqrt(2d0*dtau*diffcHe)*gasdev(idum)!*0d0 
       Phiz=dsqrt(2d0*dtau*diffcHe)*gasdev(idum)!*0d0 

!       !1a-Se rota el helio con respecto al ejex de laboratorio un angulo phix 
       call rotacion(ei(1),ei(2),ei(3),rvec_he(1),rvec_he(2),rvec_he(3),phix) 
!        
!       !1b-Se rota el helio con respecto al ejey de laboratorio un angulo phiy 
       call rotacion(ej(1),ej(2),ej(3),rvec_he(1),rvec_he(2),rvec_he(3),phiy) 
!        
!       !1c-Se rota el helio con respecto al ejez de laboratorio un angulo phiz 
       call rotacion(ek(1),ek(2),ek(3),rvec_he(1),rvec_he(2),rvec_he(3),phiz) 

      !if (abs(rdis(xcor(:,j),zcor(:,j,jhe))-rfb).gt.tol) stop 'distancia'
             enddo

        else        

!       !primero, se mueve el centro de masa de la molecula y el atomo de helio
       do i=1,3 
       xcor(i,j)=xcor(i,j)+gasdev(idum)*dsqrt(dtau*2d0*diffccm)+Fq_x(i) 
       do jhe=1,nhe
       zcor(i,j,jhe)=zcor(i,j,jhe)+gasdev(idum)*dsqrt(dtau*2d0*diffcHe)+Fq_z(i,jhe) 
              ! write(*,*) 'Atilahe',dtau,diffche
              ! stop
       enddo
       enddo
       endif

!       !Ahora se calculan los angulos para rotar la molecula
       Phix=dsqrt(2d0*dtau*A)*gasdev(idum)!*0d0 
       Phiy=dsqrt(2d0*dtau*B)*gasdev(idum)!*0d0 
       if (moltyp.gt.2) Phiz=dsqrt(2d0*dtau*C)*gasdev(idum)!*0d0 
       
!       !Se llama a la subrutina que efectua la rotacion 
!        
!       !2a-Se rota la molecula con respecto al ejex molecular un angulo phix 
       call rotacion(y1(1,j),Y1(2,j),Y1(3,j),y3(1,j),y3(2,j),y3(3,j),phix) 
       call rotacion(y1(1,j),Y1(2,j),Y1(3,j),y2(1,j),y2(2,j),y2(3,j),phix) 
!        
!       !2b-se rota la molecula con respecto al ejey molecular un angulo phiy 
       call rotacion(y2(1,j),Y2(2,j),Y2(3,j),y3(1,j),y3(2,j),y3(3,j),phiy) 
       call rotacion(y2(1,j),Y2(2,j),Y2(3,j),y1(1,j),y1(2,j),y1(3,j),phiy) 
!        
!       !2c-se rota la molecula con respecto al ejez molecular un angulo phiz 
       if (moltyp.gt.2) then
       call rotacion(y3(1,j),Y3(2,j),Y3(3,j),y1(1,j),y1(2,j),y1(3,j),phiz) 
       call rotacion(y3(1,j),Y3(2,j),Y3(3,j),y2(1,j),y2(2,j),y2(3,j),phiz)
       endif
      
   

       enddo


       return
       end Subroutine diffusion

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Subrutina que calcula el potencial de cada una de las réplicas!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Subroutine potential_calculation(npaso,N,icontr)
        real (rk) :: phi,theta,psi,rhehe,vvhehe,v_molhe,sumpot      
        real (rk) :: phig,thetag,rg,xvec(3)      
        real (rk),dimension(3)::rrot,rrotb
        real(rk):: rhe(3,nw+delta_N_max,nhe)
        integer(ik) :: j,jhe,khe,icontr,npaso
        integer(ik),intent(in)::N 
        real(rk)::xcmass(3)


        do j=1,N
        !Euler angles calculation
      ! call euler(y1(:,j),y2(:,j),y3(:,j),phi,theta,psi)


        sumpot=0.0d0
          do jhe=1,nhe
            !rhe(:,j,jhe)=(zcor(:,j,jhe)-xcor(:,j))*bo2ar
             !Molecule is aligned with lab axis
            !call rotaciones(rhe(:,j,jhe),rrot,phi,theta,psi)
             call potcalc(y1(:,j),y2(:,j),y3(:,j),xcor(:,j),zcor(:,j,jhe),v_molhe)
             sumpot=sumpot+v_molhe/har2cm
          enddo

          pot(j)=sumpot

          sumpot=0d0

       !A este potencial hay que sumarle todas las interacciones he-he

       do jhe=1,nhe
        do khe=jhe+1,nhe

          rhehe=dsqrt(dot_product(zcor(:,j,jhe)-zcor(:,j,khe),zcor(:,j,jhe)-zcor(:,j,khe)))

          call pothehe(vvhehe,rhehe)
        
          sumpot=sumpot+vvhehe

        enddo
       enddo

      pot(j)=pot(j)+sumpot

      if (state.ne.0) then
      call wvfunct(phi,theta,psi,j,mod(npaso,2))
      endif
      signo(j)=fonda(j,1)*fonda(j,0)

       if (npaso.ge.nstps-100/nruns) call dist_writing(npaso,rrot)
       Enddo 
      
        return 
        End Subroutine potential_calculation 



        subroutine dist_writing(npaso,rrot)
          integer(ik),intent(in)::npaso      
          real(rk),intent(in)::rrot(3)
          real(rk)::rg,thetag,phig 
      !Creacion del archivo de distribuciones
       rg=dsqrt(dot_product(rrot,rrot))
       thetag=dacos(rrot(3)/rg)
       phig=datan(rrot(2)/rrot(1))
       write(250,1000) rrot(1),rrot(2),rrot(3),rg,phig,thetag

1000 format(7(F10.4))        





                 return
                end Subroutine dist_writing
         
         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !Subrutina que calcula las particulas que mueren, las que  
        !sobreviven y las que nacen 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        Subroutine branch(npaso,N,icorrida) 
        real (rk) :: xnew(3,nw+delta_N_max),y1new(3,nw+delta_N_max),y2new(3,nw+delta_N_max)
        real (rk) ::y3new(3,nw+delta_N_max),znew(3,nw+delta_N_max,nhe),potnew(nw+delta_N_max)
        integer :: npaso,N,mbranch,nn,ikount,j,jhe,icorrida,kk 
        real (rk) :: W,alphar,erme_wn2,er_wn,erme,term,ekineticmol,ekinetiche,erme2
        real (rk) :: sumpotential,fondanew(nw+delta_N_max),aleat

       ikount=0 
!       !Este ciclo calcula mn para cada particula y suma el potencial total 
!       !para todos los caminantes vivos 
!         
       do j=1,N 

       If (is) then
       Call Ecinetica(j)
       Else
       Eloc(j)=0d0
       Endif

       Eloc(j)=Eloc(j)+pot(j) 

       W=dexp(-(Eloc(j)-Er(npaso-1))*dtau) 
!        
       aleat=ran2(idum)
       If(W+aleat.gt.3d0) Then 
       mbranch=3 
       Else 
       mbranch=Int(W+aleat) 
       Endif 

       If (signo(j).lt.0d0) mbranch=0 

!       If (mod(npaso,2).ne.0) then
!       nnn=1
!       Else
!       nnn=2
!       endif


        
                                 
             ikount=ikount+mbranch

          !Si se excede el numero maximo de caminantes permitido, se para 
            If (ikount.gt.nw+delta_N_max) Then 
            stop 'Mas caminantes que los permitidos'  
            endif

             do nn=ikount-mbranch+1,ikount  
             potnew(nn)=pot(j) 

            kk=mod(npaso,2)
            fondanew(nn)=fonda(j,kk) 

            xnew(:,nn) = xcor(:,j)  
           

           do jhe=1,nhe  
            znew(:,nn,jhe) = zcor(:,j,jhe)  
           enddo
                                               
            y1new(:,nn) = y1(:,j)  
            y2new(:,nn) = y2(:,j)  
            y3new(:,nn) = y3(:,j)  
            enddo
                        
             Enddo 

!        !Este ciclo asigna los valores nuevos de todos las variables a los nuevos caminantes 
        do j=1,ikount 
!            
            pot(j)=potnew(j)
            fonda(j,kk)=fondanew(j)
            
            xcor(:,j) = xnew(:,j)  
           

          do jhe=1,nhe	  
            zcor(:,j,jhe) = znew(:,j,jhe)  
          enddo  
                                                  
            y1(:,j) = y1new(:,j)  
            y2(:,j) = y2new(:,j)  
            y3(:,j) = y3new(:,j)  

         enddo 
         
         
         
        !Aqui se calcula la nueva energia de referencia para el proximo paso 
          
          sumpotential=0d0 

          do j=1,ikount 
          sumpotential=sumpotential+pot(j) 
          enddo 

          sumpotential=sumpotential/dfloat(ikount) 

          N=ikount 
        
         term = dfloat(N)/dfloat(nw) 
         alphar = 1.d0/dtau 
         er(npaso) = sumpotential + alphar*(1.d0-term)         
         ermedia2=ermedia2+er(npaso)
         erme2=ermedia2/dfloat(npaso)
         erme_wn2=erme2*har2cm 


         if (npaso.gt.dfloat(nstps/2)) Then 
         ermedia=ermedia+er(npaso) 
         erme=ermedia/(dfloat(npaso)-dfloat(nstps/2)) 
         erme_wn(icorrida)=erme*har2cm 
         endif 
         er_wn=er(npaso)*har2cm 
         if (Abs(dfloat(Int(npaso/100))-dfloat(npaso)/100d0).lt.1d-10) Then 
         write(6,55) npaso,er_wn,erme_wn(icorrida),erme_wn2,ikount,icorrida 
         endif 

       55 format (I6,2x,3(F14.8,2x),I6,2x,I2)
        
          
        return 
        end Subroutine branch 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Subroutine write_emedia
   real (rk):: dato,media,sd
   integer:: j
   
       media=0d0
       sd=0d0
       
       do j=1,nruns
       media=media+erme_wn(j)
       enddo
       media=media/nruns
       do j=1,nruns
       sd=sd+(erme_wn(j)-media)**2
       enddo       
       sd=sqrt(sd/(nruns-1))
       
       write (300,1000) media,sd,nruns

1000  Format(2(F14.6,2x),I2)

       return 

       End subroutine write_emedia
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
        subroutine wvfunct(phi,theta,psi,j,ind)
        real (rk)::phi,theta,psi
        integer::ind,j
         
         !fonda(j,ind)=0d0
         if (state.eq.1) then
         fonda(j,ind)=dcos(theta)
         return
         endif
         
         if (state.eq.2) then
         fonda(j,ind)=dsin(theta)*dcos(phi)
         return
         endif
     
        end subroutine wvfunct
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        subroutine rotacion(xn,yn,zn,xh,yh,zh,phiang) 
        real (rk)::r,c,s,t,x,y,z,xn,yn,zn,r1,r2,r3,xh,yh,zh,phiang
        
        !  rotation from  
        !c     ` 
         
              r=dsqrt(xn**2+yn**2+zn**2) 
              
              c=dcos(phiang) 
              s=dsin(phiang) 
              t=1.d0-c 
               
        !c     unit vector along axis of rotation 
         
              x=xn/r 
              y=yn/r 
              z=zn/r 
               
        !c     point to be rotated 
         
              r1=xh 
              r2=yh 
              r3=zh 
         
              xh=c*r1+s*(z*r2-y*r3)+x*t*(x*r1+y*r2+z*r3) 
              yh=c*r2+s*(x*r3-z*r1)+y*t*(x*r1+y*r2+z*r3) 
              zh=c*r3+s*(y*r1-x*r2)+z*t*(x*r1+y*r2+z*r3) 
         
              return 
              end Subroutine rotacion 
         
         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Esta subrutina rota un vector en el labframe a un vector en el sistema rotado usando los! 
        !ángulos de euler.                                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Subroutine  rotaciones(y1,y2,phi,theta,psi)
          Real (rk),intent(in):: y1(3),phi,theta,psi
          Real (rk),intent(inout)::y2(3)
          Real (rk)::matb(3,3)
          Real (rk)::prueba1,prueba2
          
              !Comprobación. Se realiza la transformación para alinear los vectores
                  
                
                      !Matriz directa, al multiplicar las coordenadas de un vector en el sistema
                      !fijo, nos da las coordenadas del mismo en el sistema rotado.

                     matb=reshape((/cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi),&
                                   -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi),&
                                    sin(theta)*sin(phi),&
                                    cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi),&
                                   -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi),&
                                   -sin(theta)*cos(phi),&
                                    sin(theta)*sin(psi),&
                                    sin(theta)*cos(psi),&
                                    cos(theta)/),(/3,3/))
                    
               y2=matmul(matb,y1)
         
             
            return
            end subroutine rotaciones


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Subrutina que calcula el potencial he-he
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              subroutine pothehe(vv1,rhe) !he-he potential
              Real (rk) r,vv1,rm,rhe,x,vst1,vst2,fx,vst,econ,alpha,beta,aa,d,eps,c6,c8,c10
             
              aa=1.8443101d05
              alpha=10.43329537
              beta=-2.27965105
              rm=2.963
              d=1.4826
              eps=10.948
              econ=3.157865d05
              c6=1.36745214
              c8=0.42123807
              c10=0.17473318
              r=0.529177*rhe
              x=r/rm
              vst1=aa*dexp(-alpha*x+beta*x**2)
              vst2=c6/x**6+c8/x**8+c10/x**10
              if(x.lt.d) then
                 fx=dexp(-(d/x-1.d0)**2)
              else
                 fx=1.d0
              endif
              vst=vst1-fx*vst2
              vv1=eps*vst/econ
              return
              end subroutine pothehe 
              

           

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !Subrutina que calcula la fuerza cuantica para molecula y atomos de helio!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE Fq (j,Fq_z,Fq_x)
               INTEGER :: j,k, kk, kkk
               REAL (rk) ::  R
               REAL (rk),  DIMENSION (3,nhe) :: FDER
               REAL (rk),  DIMENSION (3,nhe) :: Fq_z
               REAL (rk),  DIMENSION (3) :: Fq_x
               REAL (rk),  DIMENSION (3) :: Rcords
              

               DO k=1,nhe
                        Rcords=xcor(:,j)-zcor(:,j,k)
                        R=DSQRT(DOT_PRODUCT(Rcords,Rcords))

                        DO kk=1,3 
                        FDER(kk,k)=(xcor(kk,j)-zcor(kk,j,k))*(5.d0*cc/R**7-aa/R)
                        ENDDO
               ENDDO

               DO k=1,nhe
                  DO kk=1,3
                  Fq_z(kk,k)=-dtau/mHe*FDER(kk,k)
                  ENDDO   
               ENDDO

              
                  DO kk=1,3
                  Fq_x(kk)=dtau/mtot*SUM(FDER(kk,:))
                  ENDDO   
               
               return
                END Subroutine  Fq
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !subrutina que calcula la parte de la energia cinetica para elocal. Posteriormente se suma el termino del potencial en      
        !el programa que llama a esta subrutina.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Subroutine Ecinetica(j)
        real(rk) ::  R,kinetic,Rkkk,suma,ekineticmol,ekinetiche
        real(rk) ::  dd1,dd2
        real(rk),DIMENSION (3) :: Rcords,rkk
        real(rk),Dimension(nhe)::d2d
        integer::k,kk,kkk,j

        ekineticmol=0d0
        ekinetiche=0d0

       do k=1,nhe
           Rcords=xcor(:,j)-zcor(:,j,k) 
           R=DSQRT(DOT_PRODUCT(Rcords,Rcords))  
            
       do kk=1,3
        suma=0d0
        

        do kkk=1,nhe
            rkk=xcor(:,j)-zcor(:,j,kkk) 
            Rkkk=DSQRT(DOT_PRODUCT(rkk,rkk)) 
           suma=suma+dfx(rkkk,kk,kkk,j)
        enddo  

        dd2=d2fx(r,kk,k,j)
        dd1=dfx(r,kk,k,j)
        ekineticmol=ekineticmol+dd2+dd1*suma
        ekinetiche=ekinetiche+dd2+dd1**2 


       enddo
       enddo
          

       Eloc(j)=-1d0/(2d0*mHe)*ekinetiche-1d0/(2d0*mtot)*ekineticmol


        return 
        end subroutine Ecinetica
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function df(rr)
        Real (rk)::rr

        df=5d0*cc/rr**6-aa

        endfunction df
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function d2f(rr)
        Real (rk)::rr

        d2f=-30d0*cc/rr**7

        endfunction d2f


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function dfx(rr,kder,jhe,j)
        Real (rk)::rr
        integer::kder,jhe,j

        dfx=df(rr)*(xcor(kder,j)-zcor(kder,j,jhe))/rr

        endfunction dfx

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Real (rk) Function d2fx(rr,kder,jhe,j)
           Real (rk)::rr,dif
          integer::kder,jhe,j

        dif=xcor(kder,j)-zcor(kder,j,jhe)
        d2fx=d2f(rr)*(dif**2)/rr**2+df(rr)*(1-dif**2/rr**2)/rr

        endfunction d2fx


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine euler(eir,ejr,ekr,phi,theta,psi)
         Real (rk):: matb(3,3)
         Real (rk):: eirr(3),ejrr(3),ekrr(3),cv(3)
         Real (rk), Intent(in)::eir(3),ejr(3),ekr(3)
         Real (rk),parameter::ei(3)=(/1d0,0d0,0d0/),ej(3)=(/0d0,1d0,0d0/),ek(3)=(/0d0,0d0,1d0/) 
         Real (rk)::phi,theta,psi,pru1,pru2,pru3,pru4,pru5,pru6,pi
         Real (rk)::phip,thetap,psip,valprueba1,valprueba2  
         Real (rk)::ter11,ter12,ter13,ter21,ter22,ter23,ter31,ter32,ter33,prudef 
         integer::ind1,ind2 


            pi=dacos(-1d0)

           
                !Vectores rotados           
                eirr=eir
                ejrr=ejr
                ekrr=ekr

                 !Calculo de los angulos de Euler

               
         
                 theta=dacos(ekrr(3))
                 if (theta.gt.1d-14) then
                 phi=datan2(ekrr(1),-ekrr(2))
                 if (phi.lt.0d0) phi=phi+2d0*pi
                 !print*,  'componentes Euler',eirr(3),e
                 psi=datan2(eirr(3),ejrr(3)) 
                 if (psi.lt.0d0) psi=psi+2d0*pi
                 else
                 phi=dacos(eirr(1))
                 psi=0d0
                 endif
              
            
                !Comprobación. Se realiza la transformación para alinear los vectores
                  
                
                      !Matriz directa, al multiplicar las coordenadas de un vector en el sistema
                      !fijo, nos da las coordenadas del mismo en el sistema rotado.

                     matb=reshape((/cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi),&
                                   -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi),&
                                    sin(theta)*sin(phi),&
                                    cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi),&
                                   -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi),&
                                   -sin(theta)*cos(phi),&
                                    sin(theta)*sin(psi),&
                                    sin(theta)*cos(psi),&
                                    cos(theta)/),(/3,3/))

                    
                   !Inversa de la matriz anterior, al multiplicar las coordenadas de un vector en el sistema
                   !rotado, nos da las coordenadas del mismo en el sistema fijo.



                  eirr=Matmul(Matb,eirr)
                  ejrr=Matmul(matb,ejrr)
                  ekrr=Matmul(matb,ekrr)

             
              
                 pru1=dot_product(ei-eirr,ei-eirr)
                 pru2=dot_product(ej-ejrr,ej-ejrr)
                 pru3=dot_product(ek-ekrr,ek-ekrr)
                 
                 prudef=pru1+pru2+pru3

                if (prudef.gt.1d-15) then
                print*,'mal algo en angulos de Euler'
                print*,'pruebas',pru1,pru2,pru3,prudef
                print*,'angulos',phi,theta,psi
                stop
                endif
                     

                  end Subroutine euler
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            function rdis(a,b)
            real(rk)::rdis,a(3),b(3),c(3)

            c=b-a
            rdis=dsqrt(dot_product(c,c))
               
            return
            end function rdis

            end Module dmc_cal
