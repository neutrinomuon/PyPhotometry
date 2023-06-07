***** D. Le Borgne : this file is identical to the PEGASE.2 's calib.f file
*     Except that the path where the filters.dat fit can be found is set 
*     by the variable PEG_ROOT.

***** Calibration of the filters

      program calib_HR

      implicit none

      !include 'peg_config.f'    ! to find where "data" are
      !include 'peg_include.f' 

      integer nmaxfilters,nmaxlambdatrans
      parameter (nmaxfilters=10000,nmaxlambdatrans=10000)

      character*72 name(nmaxfilters),arq_lix
      integer nlambdaVega,nfilters,nlambdafilter(nmaxfilters),i,j
      integer typetrans(nmaxfilters),nlambdaFSD,nlambdaSun
      integer typecalib(nmaxfilters)
      double precision lambdafilter(nmaxfilters,nmaxlambdatrans)
      double precision trans(nmaxfilters,nmaxlambdatrans)
      double precision lambdamed,transmed,fluxtrans
      double precision area(nmaxfilters),lambdaeff,fluxtransFSD
      double precision fluxVega(5000),lambdaVega(5000)
      double precision lambdamean(nmaxfilters)
      double precision ABVega,a,TGVega,lambdaeffFSD
      double precision areanu(nmaxfilters),lambdaFSD(5000)
      double precision fluxFSD(5000)
      double precision lambdaSun(5000),fluxSun(5000)
      double precision lambdaeffSun,fluxtransSun,c
      parameter(c=2.99792458d18)
      

***** Reading of the spectrum of Vega

      open(10,status='old',file='VegaLR.dat')
      read(10,*) arq_lix,nlambdaVega
      do i=1,nlambdaVega
         read(10,*) lambdaVega(i),fluxVega(i) 
      end do
      close(10)

***** Reading of the spectrum of the Sun

      open(10,status='old',file='Sun_LR.dat') 
*     Intrinsic flux: erg/s/A
      read(10,*) arq_lix,nlambdaSun     
      do i=1,nlambdaSun
         read(10,*) lambdaSun(i),fluxSun(i) 
      end do
      close(10)

***** Reading of the spectrum of BD+17o4708 (F subdwarf used to calibrate the 
***** Thuan & Gunn system).

      open(10,status='old',file='BD+17o4708.dat')
      read(10,*) arq_lix,nlambdaFSD     
      do i=1,nlambdaFSD
         read(10,*) lambdaFSD(i),fluxFSD(i) 
      end do
      close(10)

***** Reading of the transmissions of the filters

      open(20,status='old',file='filters.dat')
      read(20,*) nfilters

      !nfilters = 1
          
      do i=1,nfilters
         read(20,*) nlambdafilter(i),typetrans(i),typecalib(i),name(i)
c        write (*,'(a,3(i15),a)') '#',nlambdafilter(i),typetrans(i)
c    $        ,typecalib(i),name(i)
         do j=1,nlambdafilter(i)
            read(20,*) lambdafilter(i,j),trans(i,j)
c            write (*,*) lambdafilter(i,j),trans(i,j)
         end do        
      end do
      close(20)

      open(40,status='unknown',file='calib.dat')

      write(40,'(a,4x,a,4x,a,3x,a,4x,a,2x,a,2x,a,2x,a,2x,a)') 
     $     'filter','index','Flambda(Vega)',
     $     'area','mean lambda','lambda eff(Vega)','mAB(Vega)',
     $     'mTG(Vega)','Flambda(Sun)'

***** Areas and mean wavelengths of the filters

      do i=1,nfilters
         area(i)=0.
         a=0.
         areanu(i)=0.
         fluxtrans=0.
         lambdamean(i)=0.
         do j=1,nlambdafilter(i)-1
            lambdamed=(lambdafilter(i,j)+lambdafilter(i,j+1))/2. 
            transmed=(trans(i,j)+trans(i,j+1))/2.
            lambdamean(i)=lambdamean(i)+lambdamed**(1+typetrans(i))
     $           *transmed*(lambdafilter(i,j+1)-lambdafilter(i,j))
            area(i)=area(i)+transmed
     $           *(lambdafilter(i,j+1)-lambdafilter(i,j))
     $           *lambdamed**typetrans(i)
            areanu(i)=areanu(i)+transmed
     $           *(lambdafilter(i,j+1)-lambdafilter(i,j))
     $           *lambdamed**typetrans(i)/lambdamed**2*c

c            write (*,*) lambdamed,transmed,lambdamean(i),area(i)
c     $           ,areanu(i)


         end do
         lambdamean(i)=lambdamean(i)/area(i)
 
         write (*,'(i5,2x,a25,2(e15.5))') typetrans(i),name(i)
     $        ,lambdamean(i),areanu(i)
   
      end do

***** Computation of the effective wavelengths, AB-magnitudes of Vega
***** and calibration constants of the filters.

      do i=1,nfilters
         call calculflux(lambdaVega,fluxVega,lambdafilter,trans,area,
     $        nlambdafilter,typetrans,fluxtrans,lambdaeff,i,nlambdaVega)
         call calculflux(lambdaSun,fluxSun,lambdafilter,trans,area,
     $        nlambdafilter,typetrans,fluxtransSun,lambdaeffSun,i,
     $        nlambdaSun)
         ABVega=-2.5*log10(fluxtrans*area(i)/areanu(i))-48.60
         call calculflux(lambdaFSD,fluxFSD,lambdafilter,trans,area,
     $        nlambdafilter,typetrans,fluxtransFSD,lambdaeffFSD,i,
     $        nlambdaFSD)
         if (fluxtransFSD.gt.0.) then
            TGVega=-2.5*log10(fluxtrans/fluxtransFSD)+9.5
         else
            TGVega=99.999
         end if
         write(40,50) name(i),i,fluxtrans,
     $        area(i),lambdamean(i),lambdaeff,
     $        ABVega,TGVega,fluxtransSun
      end do
      close(40)

 50   format(a10,1x,i3,7x,e9.4,3x,e9.4,3x,e9.4,5x,e9.4,5x,
     $     2(f8.3,3x),2x,e9.4)

      end

*****

      subroutine calculflux(lambda,flux,lambdafilter,trans,area,
     $     nlambdafilter,typetrans,fluxtrans,lambdaeff,i,nlambda)

      implicit none
      !include 'peg_include.f'

      integer nmaxfilters
      parameter (nmaxfilters=10000)

      integer nlambda
      integer i,j,k,nlambdafilter(nmaxfilters),typetrans(nmaxfilters)
      double precision fluxtrans,lambda(5000),flux(5000),transsup
      double precision trans(nmaxfilters,1000),fluxinf,transinf,fluxsup
      double precision lambdaeff,lambdainf,lambdasup
      double precision area(nmaxfilters),lambdafilter(nmaxfilters,1000)
      double precision interplinlin,interploglog

      lambdaeff=0.
      fluxtrans=0.
      if ((lambda(1).gt.lambdafilter(i,1)).or.(lambda(nlambda).lt.
     $     lambdafilter(i,nlambdafilter(i)))) then
         fluxtrans=-1.
      else
         k=1
         do while(lambda(k+1).lt.lambdafilter(i,1))
            k=k+1
         enddo
         j=1
         lambdainf=lambdafilter(i,j)
         transinf=0.
         fluxinf=0.
         lambdasup=lambdainf
         do while(j+1.le.nlambdafilter(i).and.k+1.le.nlambda)
            if (lambdafilter(i,j+1).lt.lambda(k+1)) then
               lambdasup=lambdafilter(i,j+1)

               !write (*,*) 'É agora'

               fluxsup=interploglog(lambda(k),lambda(k+1),flux(k),
     $              flux(k+1),lambdasup)
               fluxsup=interplinlin(lambda(k),lambda(k+1),flux(k),
     $              flux(k+1),lambdasup)          
               transsup=trans(i,j+1)*lambdasup**typetrans(i)               

               !stop
               j=j+1
            else
               lambdasup=lambda(k+1)
               fluxsup=flux(k+1)
               transsup=interplinlin(lambdafilter(i,j),
     $              lambdafilter(i,j+1),trans(i,j),trans(i,j+1),
     $              lambdasup)*lambdasup**typetrans(i)  
               k=k+1
            endif
            fluxtrans=fluxtrans+(transinf*fluxinf
     $           +transsup*fluxsup)*(lambdasup-lambdainf)/2
            lambdaeff=lambdaeff+(transinf*fluxinf
     $           +transsup*fluxsup)*(lambdasup-lambdainf)/2
     $           *(lambdasup+lambdainf)/2
            lambdainf=lambdasup
            fluxinf=fluxsup
            transinf=transsup
         enddo
         lambdaeff=lambdaeff/fluxtrans
         fluxtrans=fluxtrans/area(i)             
      endif

      end

***** 
      double precision function interplinlin(x1,x2,y1,y2,t)

      implicit none
      double precision x1,x2,y1,y2,t

      interplinlin=y1+(y2-y1)*(t-x1)/(x2-x1)

      !write (*,*) interplinlin

      end

***** 

      double precision function interploglog(x1,x2,y1,y2,t)

      implicit none
      double precision x1,x2,y1,y2,t,eps
      parameter(eps=1.d-37)

      if (y1.gt.eps.and.y2.gt.eps) then      
         interploglog=y1*(y2/y1)**(log(t/x1)/log(x2/x1))
      else
         interploglog=0.d0
      endif

      !write (*,*) interploglog

      end













