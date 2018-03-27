Module hcldimpot
use types
implicit none
real(rk),parameter:: pi=dacos(-1.0_rk)
real(rk):: lnfac(0:1000)
integer(ik):: evenodd(-2000:2000)
real(rk) facto(0:100)
real(rk):: rhcl= 2.415069d0     
real(rk):: a   = 3.07d0         
real(rk):: b   = 186d0          
real(rk):: c1  = 2.412008d0     
real(rk):: c2  = 7.46886431d0   
real(rk):: c3  = 1.0d0          
real(rk):: ksr(2)=(/ 0.0170d0,0.00979d0/)      
real(rk):: q(0:5)= (/0.0_rk,1.08716_rk,3.48882533_rk,2.03363010_rk,0.67488208_rk,4.26001465_rk/)
real(rk):: k(0:11,0:5,0:5,0:7)=0.0_rk
integer(ik)::isubv(31)=(/0,1,2,6,11,1,3,7,9,1,6,3,9,4,7,4,10,1,9,5,9,9,9,7,10,1,1,7,10,7,8/) 
integer(ik)::lav(31)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3/) 
integer(ik)::lbv(31)=(/0,0,0,0,0,1,1,1,1,2,2,3,3,4,5,1,2,2,2,2,2,3,3,3,3,4,4,4,4,5,3/)
integer(ik)::lv (31)=(/0,0,0,0,0,1,1,1,1,2,2,3,3,4,5,0,1,3,3,2,4,1,3,5,5,2,4,6,6,7,6/)
real(rk)::kass(31)=(/0.0_rk,0.02270071_rk,0.02136_rk,-10000._rk,0.0_rk,-0.00525001_rk,&
                    &-0.000498_rk,11218._rk,0.0_rk,-0.00221_rk,-543._rk,-0.00048769_rk,&
                    &0.0_rk,0.0000300_rk,191._rk,-0.0000499_rk,0.0_rk,-0.000241_rk,&
                    &0.0_rk,0.0000084_rk,0.0_rk,0.0_rk,0.0_rk,286._rk,0.0_rk,0.000137_rk,&
                    &0.000122_rk,343._rk,0.0_rk,-243._rk,1459._rk/)
real(rk),parameter::htocm=219474.6354_rk,cf=43.597482_rk,bohtoar=0.529177249_rk
 


Contains

       subroutine potkass
        integer(ik)::i
         
         do i=1,31
           k(isubv(i),lav(i),lbv(i),lv(i))=kass(i)
         enddo   

!     calculate factorials for apleg routine
       facto(0)=1d0
       do 5 i=1,30
 5      facto(i)=fact(i)

!     set up factorials (for 3,6,9-j symbols) and parity function
         lnfac(0) = 0.d0
         lnfac(1) = 0.d0
         do 20 i=2,1000
            lnfac(i) = lnfac(i-1) + dlog(dfloat(i))
 20      continue
         do 30 i=0,1000
            evenodd(i) = (-1)**i
            evenodd(-i) = evenodd(i)
 30      continue
       return
       end subroutine potkass

      subroutine pothcl2(rr,costh1,costh2,taud,totalcm)
!     generate potential 
!     Matt Elrod March 1993
      implicit none
      real(rk),intent(in)::rr,costh1,costh2,taud
      real(rk)::totalcm
      integer jhclmax,potqn(200,3),nterms
      real(rk) bhcl,m1,m2,rm,r
      real(rk) av(1000)
      real(rk) th1,th2,prefac,zerom,tzerom,rau,tan1p
      real(rk) tan1m,tan2m,total,total1,what1,what2,y1,y2,y3 
      real(rk) tau,totpot,this,d1,d2,d3,elec,elecnot
      real(rk) re,gamma,qpt(64),qwt(64)
      integer(ik)::i,icount,isub,la,lb,lm,m,l
      save

         r=rr
         rau=rr/bohtoar
        !th1=th1d*pi/180d0
        !costh1=dcos(th1)
        !th2=th2d*pi/180d0
        !costh2=dcos(th2)
        !tau=taud*pi/180d0

         d1=rhcl-c1
         d2=rhcl-c1
         d3=rau-c2
         y1=1d0 - dexp(-c3*d1)
         y2=1d0 - dexp(-c3*d2)
         y3=1d0 - dexp(-d3)
         tan1p=0.5d0*(1d0+dtanh(2.d0*(rau-6.d0)))
         tan1m=0.5d0*(1d0-dtanh(2.d0*(rau-6.d0)))
         tan2m=0.5d0*(1d0-dtanh(0.5d0*(rau-8.5d0)))

!     calculate expansion coefficients
         icount=1
         totpot=0d0
!     A(0,0,0) term
         av(1)=k(0,0,0,0) + (ksr(1)*dexp(-d3)&
              & +ksr(2)*dexp(-2d0*d3))*tan1m + (k(1,0,0,0)*dexp(-d3)&
              & +k(2,0,0,0)*dexp(-2d0*d3))*tan1p*tan2m&
              & +(k(6,0,0,0)*tan1p)/(rau**6d0)& 
              & +k(11,0,0,0)*(y1**2d0+y2**2d0)

!     anisotropic A terms
         do 60 la=0,5
            do 70 lb=la,5
               if ((la.eq.lb).and.(la.eq.0)) goto 70
               do 80 l=abs(la-lb),(la+lb)
!     nonelectrostatic terms
                 elecnotcal: select case(l) 
                 case(8:)
                 elecnot=0d0
                 case(0:7)
                  elecnot=(k(1,la,lb,l)*dexp(-d3)&
                  &+k(2,la,lb,l)*dexp(-2d0*d3) + k(3,la,lb,l)*dexp(-3d0*d3)&
                  &+k(4,la,lb,l)*dexp(-4d0*d3) + k(5,la,lb,l)*dexp(-5d0*d3)&
                  &+k(9,la,lb,l)*y1 + k(10,la,lb,l)*y2)*tan1p*tan2m&
                  &+(k(6,la,lb,l)/(rau**6d0) + k(7,la,lb,l)/(rau**7d0)&
                  &+k(8,la,lb,l)/(rau**8d0))*tan1p
                 case default
                    stop 'angular momentum l out of range'
                  end select elecnotcal 

!     plus electostatic terms
!     cf is a conversion factor to go from Debye and Angstroms to hartrees
                  elec=0d0
                  if (((la+lb).eq.l).and.(l.ne.0)) then
!     tan1p=1d0
              elec = evenodd(lb)*8d0*(pi**1.5d0)*((1d0/(2*l+1))**0.5d0)&
          &*(dexp(lnfac(2*la+2*lb)-lnfac(2*la+1)-lnfac(2*lb+1)))**0.5d0&
                      &*q(la)*q(lb)*tan1p/(cf*r**(l+1))
                  endif
!     remember to change tan1p
                  total=elecnot+elec


!     coefficient must be greater than 1 cm-1
                  if (abs(total*htocm).ge.1d-6) then
                     icount=icount+1
                     potqn(icount,1)=la
                     potqn(icount,2)=lb
                     potqn(icount,3)=l
                     av(icount)=total


!     A(la,lb,l)=-A(lb,la,l)
                     if (la.ne.lb) then
                        icount=icount+1
                        potqn(icount,1)=lb
                        potqn(icount,2)=la
                        potqn(icount,3)=l
                        av(icount)=evenodd(la+lb)*av(icount-1)
                     endif
                  endif
                  
 80            continue
 70         continue
 60      continue


         nterms=icount
         total=0d0
         total1=0d0
!     sum over potential basis functions
         do 100 i=1,nterms
            la=potqn(i,1)
            lb=potqn(i,2)
            l=potqn(i,3)
            prefac=evenodd(la-lb)*((2d0*l+1d0)/(2d0*pi**0.5d0))
            zerom=threej(la,lb,l,0,0,0)*apleg(la,0,costh1)&
                 &*apleg(lb,0,costh2)
            tzerom=0d0
            lm=la
            if (la.gt.lb) lm=lb 
            do 85 m=1,lm
               tzerom=tzerom + evenodd(m)*2d0*threej(la,lb,l,m,-m,0)&
                   &*apleg(la,m,costh1)*apleg(lb,m,costh2)*cos(m*tau)
               what1=apleg(la,m,costh1)
               what2=apleg(lb,m,costh2)
 85         continue
          total=total+prefac*(zerom+tzerom)*av(i)

 100     continue
         totalcm=total*htocm
         return
         end subroutine pothcl2 

      function apleg(l,m,x)
!     this function evaluates the associated legendre polynomials for 
!     arbitrary j and k>0 or =0.  the routine follows the example from numerical 
!     recipes by Press, Flannery, Teukolsky, and Vetterling.  I''ve added
!     an additional factor so that the functions are normalized.
!     reproduces pleg(l,x) to the accuracy of the computer for k=0
      implicit none
      real(rk)::apleg,x,fct,pll,pmm,pmmp1,somx2,anorm
      integer(ik)::m,l,i,ll

!     write(*,*)l,m,x
      if((m.lt.0).or.(m.gt.l))then
         write(4,*)'bad args apleg',l,m
      endif

!     normalization
      if ((l.eq.0).and.(m.eq.0)) then
         anorm=(1d0/(4d0*pi))**0.5d0
      else
         anorm=(((2d0*l+1d0)*facto(l-m))/(4d0*pi*facto(l+m)))**.5d0
      endif

!     write(*,*)l,m,x,anorm
      if ((l.eq.2).and.(m.eq.2)) then
!     write(0,*)'l=2 m=2 norm =',anorm
      endif


      pmm=1d0
      if(m.gt.0)then
         somx2=dsqrt((1d0-x)*(1d0+x))
         fct=1d0
         do 11 i=1,m
            pmm=-pmm*fct*somx2
 11         fct=fct+2d0
         endif	

!     pmm refers to the l=m apleg, pmmp1 refers to l=m+1
         if(l.eq.m)then
            apleg=pmm
         else
            pmmp1=x*(2d0*m+1d0)*pmm
            
            if(l.eq.m+1d0)then
               apleg=pmmp1
            else
               do 12 ll=m+2,l
                  pll=(x*(2d0*ll-1d0)*pmmp1-(ll+m-1d0)*pmm)/(ll-m)
                  pmm=pmmp1
 12               pmmp1=pll
                  apleg=pll
               endif
            endif
            
            apleg=anorm*apleg

            return
            end function apleg


       function fact(ix)
       implicit none
       real(rk)::fact,afact
       integer(ik)::ix,i

      fact=0
      afact=1d0
      if(ix)599,601,602
 602  do 600 i=1,ix
 600     afact=afact*i
 601     fact=afact
 599     return
         end function fact



       function threej(j1,j2,j3,m1,m2,m3)
              implicit none
      integer(ik):: j1,j2,j3,m1,m2,m3,lownu,highnu,nu
      real(rk):: threej,lndelta,lnpart1,lnpart2,lnthrj,answer
      threej = 0.d0
      answer = 0.d0

!     if m1=m2=m3 then j1+j2+j3 must be even
      if ((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0).and.(mod((j1+j2+j3),2).eq.1)) return

!     check triangle rule
      if (m1+m2+m3 .ne. 0) return
      if (abs(m1).gt.j1 .or. abs(m2).gt.j2 .or. abs(m3).gt.j3) return
      if (j3 .gt. j1+j2 .or. j3 .lt. abs(j1-j2)) return

      
      lownu = max0(0,j2-j3-m1,j1+m2-j3)
      highnu = min0(j1-m1,j2+m2,j1+j2-j3)

      lndelta = .5d0*(lnfac(j1+j2-j3)+lnfac(j1+j3-j2)+lnfac(j2+j3-j1)-&
          &lnfac(j1+j2+j3+1))
      lnpart1 = .5d0*(dlog(dfloat(2*j3+1))+lnfac(j1+m1)+lnfac(j1-m1)+l&
          &nfac(j2+m2)+lnfac(j2-m2)+lnfac(j3+m3)+lnfac(j3-m3))
      lnpart1 = lnpart1 + lndelta

      do 20 nu=lownu,highnu
         lnpart2 = lnfac(j1-m1-nu)+lnfac(j3-j2+m1+nu)+lnfac(j2+m2-nu)+ &
          &   lnfac(j3-j1-m2+nu)+lnfac(nu)+lnfac(j1+j2-j3-nu)
         lnthrj = lnpart1 - lnpart2

         answer = answer + evenodd(nu)*dexp(lnthrj)
 20   continue
!     Now the V-C coefficient is done, so you just need to include the phase part 
!     and you''ll have the 3-j value.  The formula on p. 39 of Brink and Satchler 
!     gives the factors necessary to go from V-C coefficient to 3-j symbol.  They
!     use the Condon and Shortley phase convention.(I think)
      answer = answer/dsqrt(dfloat(2*j3+1))*evenodd(j1-j2-m3)
      threej = answer
      return
      end function threej

end Module hcldimpot

