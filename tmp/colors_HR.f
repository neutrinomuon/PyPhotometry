c AL- Modified from the original verion of PEGASE.2, which used to read ascii
c files. Here, fits outputs of spectra.f are read. An ascii file is created,
c as in the original version. The units of the columns are as in the output of 
c the original version of colors.f (although the flux outputs of spectra.f are 
c now in solar units rather than in the original erg/s).

c History: 
c PhP 2001/09/29: First version. Still some problems with the flux output units.
c	Copy of scenario properties still to be done.
c AL 2002/11/27 : * Units checked. Output units now as in original colors.f.
c       * fits_spec_wavl_read (in fits_spec_io.f) now reads nlines properly.
c       * Compile with colors.com (i.e. with fits_spec_io.f, util.f, 
c       and with the fitsio library).
c DLB 2005/01/19 : bug fix of D4000 measurement : area is set to 1 if typecalib.eq.0
c-----------------------------------------------------------------------------
c
c      Uses subroutines : file_open, fits_spec_open, fits_spec_wave_read,
c			  fits_spec_wavl_read, fits_spec_cont_r 
c			  fits_spec_line_r, fits_spec_para_r
c			  + fitsio library
c-----------------------------------------------------------------------------
      program colors_HR

      implicit none

      include 'peg_config.f'    ! to find where "data" are
      include 'peg_include.f'    ! constants nmax...

      character*280 filespectra,filecolors
      character*80 comment
      character*100 header
      integer i,j,nfilters,nlambdafilter(nmaxfilters),ntimes,nlambda
      integer nlines,time(nmaxotimes),iline(nmaxlines)
      integer typetrans(nmaxfilters),typecalib(nmaxfilters)
      double precision trans(nmaxfilters,1000),Lsol
      double precision lambdafilter(nmaxfilters,1000)
      double precision calib(nmaxfilters),area(nmaxfilters)
      double precision lambda(nmaxlambda),flux(nmaxlambda)
      double precision V(nmaxotimes),lambdaline(nmaxlines)
      double precision fluxline(nmaxlines),Lumline(nmaxotimes,nmaxlines)
      double precision fluxfilter(nmaxotimes,nmaxfilters)
      double precision sigmagas(nmaxotimes)
      double precision SFR(nmaxotimes),ZISM(nmaxotimes)
      double precision sigmastars(nmaxotimes),sigmaWD(nmaxotimes)
      double precision sigmaBHNS(nmaxotimes),Zstars(nmaxotimes)
      double precision Zbol(nmaxotimes),fluxIR(nmaxotimes)
      double precision agestars(nmaxotimes),agebol(nmaxotimes)
      double precision fluxbol(nmaxotimes),NLym(nmaxotimes)
      double precision mag(nmaxotimes,nmaxfilters)
      double precision nSNII(nmaxotimes),D4000(nmaxotimes)
      double precision Mbol(nmaxotimes)
      double precision fluxsol(nmaxfilters),EW(nmaxotimes,nmaxlines)
      double precision fluxVega,ABVega,TGVega,sigmasub(nmaxotimes)
      double precision nSNIa(nmaxotimes),Mgal(nmaxotimes)
      double precision tauV(nmaxotimes)

      integer lun
      integer Lfits
      integer istat,hdutype
      logical write_out
      integer nargs
      parameter(Lsol=3.826e33)

***** Reading of the transmission of the filters
 

      call file_open(
     $     PEG_ROOT//'data/user_defined/filters.dat',lun,istat)

      if(istat.ne.0) then
         write(*,*)'COLORS: could no open "filters.dat"'
         stop
      endif
      read(lun,*) nfilters
      do i=1,nfilters
         read(lun,*) nlambdafilter(i),typetrans(i),typecalib(i)
         do j=1,nlambdafilter(i)
            read(lun,*) lambdafilter(i,j),trans(i,j)
         end do         
      end do
      close(lun)

***** Reading of the calibrations of the filters

      call file_open(
     $ PEG_ROOT//'data/user_defined/calib.dat',lun,istat)

      if(istat.ne.0) then
         write(*,*)'COLORS: could no open "calib.dat"'
         stop
      endif
      read(lun,'(a)') header
      do i=1,nfilters
         read(lun,'(21x,e9.4,3x,e9.4,31x,2(f8.3,3x),2x,e9.4)') 
     $        fluxVega,area(i),ABVega,TGVega,fluxsol(i)
*     D4000 break
         if (typecalib(i).eq.0) area(i)=1.
*     ref=Vega
         if (typecalib(i).eq.1) calib(i)=2.5*log10(fluxVega)+0.03
*     ref=AB
         if (typecalib(i).eq.2) calib(i)=2.5*log10(fluxVega)+ABVega
*     ref=BD+17d4708 (Thuan & Gunn system)
         if (typecalib(i).eq.3) then
            if (TGVega.gt.99.) then
               calib(i)=99.999
               write(*,'(a,a,i3)') ' The Thuan & Gunn magnitude is',
     $              ' undefined for the filter number ',i
            else
               calib(i)=2.5*log10(fluxVega)+TGVega
            endif
         endif
*     ref=-21.10
         if (typecalib(i).eq.4) calib(i)=-21.10
*     ref=-21.175
         if (typecalib(i).eq.5) calib(i)=-21.175
      end do
      close(lun)

***** Reading arguments if present
      nargs=iargc()
      if (nargs.gt.2.or.nargs.eq.1) then 
         write(*,*) 'Error : too few or too many arguments.'
         write(*,*) 'Usage : > colors_HR spectrum_input.fits '
     &        //'colors_out.dat (command line mode)'
         write(*,*) '     or > colors_HR                   '//
     &        '(interactive mode)'
         stop
      endif
      if (nargs.eq.2) then 
         call getarg(1,filespectra)
         call getarg(2,filecolors)
      else
         write(*,'(a,$)') 'Input filename (spectra) ?'
         read(*,'(a)') filespectra
         
         write(*,'(a,$)') 'Output filename (colors) ?'
         read(*,'(a)') filecolors
         if (filecolors.eq.' ') filecolors='colors_'//filespectra
      endif

      open(50,status='unknown',file=filecolors,iostat=istat)
      if(istat.ne.0) then
         write(*,*)'COLORS: could not create output file'
         stop
      endif

***** Reading of the spectra
      call fits_spec_open(filespectra,Lfits,nlambda,ntimes,istat)
      if (istat.ne.0) then
         write(*,*)'COLORS: could not open ETS file'
         stop
      endif
      call fits_spec_wave_read(Lfits,nlambda,lambda,istat) 
      if (istat.ne.0) then
         write(*,*)'COLORS: could not get cnt wavelengths in ETS file'
         write(*,*) 'istat=',istat
         stop
      endif
      call fits_spec_wavl_read(Lfits,nlines,lambdaline,istat)
      if (istat.ne.0) then
         write(*,*)'COLORS: couldnt get line wavelengths in ETS file'
         stop
      endif

*     Copy scenario information.
      call ftmahd(Lfits,1,hdutype,istat)
      j=0
      write_out=.false.
      do i=1,55
	call ftgcrd(Lfits,'COMMENT',comment,istat)
        if (istat.ne.0) exit
        if (j.eq.1.and.comment(9:11).eq.'---') exit
        if (j.eq.1) write(50,'(a)') comment(9:80)
        if (j.eq.0.and.comment(9:11).eq.'---') j=1
      enddo
      write(50,'(a)')
     $ '************************************************************'


      write(50,*) ntimes


      do i=1,nlines
         if ((lambdaline(i)-lambda(1))*(lambdaline(i)
     $        -lambda(nlambda)).lt.0.) then
            j=1
            do while (lambda(j).lt.lambdaline(i))
               j=j+1
            end do
            iline(i)=j-1
         else
            iline(i)=0
         end if        
      end do

      do j=1,ntimes
         call fits_spec_cont_r(Lfits,nlambda,j,flux,istat)
         call fits_spec_line_r(Lfits,nlines,j,fluxline,istat)
         call fits_spec_para_r(Lfits,j,time(j),
     $        Mgal(j),sigmastars(j),sigmaWD(j),
     $        sigmaBHNS(j),sigmasub(j),sigmagas(j),
     $        ZISM(j),Zstars(j),Zbol(j),
     $        fluxbol(j),tauV(j),fluxIR(j),SFR(j),NLym(j),
     $        nSNII(j),nSNIa(j),agestars(j),agebol(j), istat)

*     Equivalent width
         do i=1,nlines
            Lumline(j,i)=fluxline(i)
            if (iline(i).ne.0.and.flux(iline(i)).gt.0.) then
               EW(j,i)=fluxline(i)/(flux(iline(i)) 
     $              +(lambdaline(i)-lambda(iline(i)))
     $              *(flux(iline(i)+1)-flux(iline(i)))
     $              /(lambda(iline(i)+1)-lambda(iline(i))))
            else
               EW(j,i)=0.
            end if
         end do

c Note: the fits_ version of spectra.f returns fluxbol in solar units 
         if (fluxbol(j).gt.0.) then
            Mbol(j)=4.75-2.5*log10(fluxbol(j))
         else            
            Mbol(j)=0.
         endif
c Note: calculflux returns fluxfilter in erg/s.
         call calculflux(j,nfilters,nlambdafilter,lambdafilter,
     $        trans,nlambda,lambda,flux,nlines,lambdaline,fluxline,
     $        area,fluxfilter,typetrans)

         do i=1,nfilters
            if (fluxfilter(j,i).gt.0.) then
               mag(j,i)=-2.5*log10(fluxfilter(j,i)/1.1965d40)
     $              +calib(i)
            else
               mag(j,i)=0.
            end if
         end do
         V(j)=mag(j,42)

*     4000A break (Bruzual 1983); at redshift=0 only.

         if (fluxfilter(j,46).gt.0..and.
     $        fluxfilter(j,47).gt.0.) then
            D4000(j)=fluxfilter(j,47)/
     $           fluxfilter(j,46) 
         else
            D4000(j)=0.
         end if
      end do
c      close(20)     
      write(50,51) 
     $     'time','Mgal','M*','MWD','MBHNS','Msub','Mgas','Zgas',
     $     '<Z*>mass','<Z*>Lbol'
      do j=1,ntimes
         write(50,'(i5,9(1x,e8.3))') time(j),Mgal(j),sigmastars(j),
     $        sigmaWD(j),sigmaBHNS(j),sigmasub(j),sigmagas(j),
     $        ZISM(j),Zstars(j),Zbol(j)
      end do
      write(50,52) 
     $     'time','Lbol','tauV','Ldust/Lbol','SFR','nSNII',
     $     'nSNIa','<t*>mass','<t*>Lbol'
      do j=1,ntimes
         write(50,'(i5,8(1x,e8.3))') time(j),fluxbol(j)*Lsol,tauV(j),
     $        fluxIR(j),SFR(j),nSNII(j),nSNIa(j),agestars(j),agebol(j)
      end do
      write(50,53) 
     $     'time','nLymcont','L(Ha)','W(Ha)','L(Hb)',
     $     'W(Hb)','LB/LBsol','LV/LVsol','D4000'
      do j=1,ntimes
         write(50,'(i5,8(1x,e8.3))') time(j),NLym(j)*Lsol,
     $        Lumline(j,2)*Lsol, EW(j,2),Lumline(j,1)*Lsol,EW(j,1),
     $        fluxfilter(j,45)/fluxsol(45),
     $        fluxfilter(j,42)/fluxsol(42),D4000(j)
      enddo
      write(50,54)
     $     'time','Mbol','V','U-B','B-V',
     $     'V-K','V-RC',
     $     'V-IC','J-H','H-K'
      do j=1,ntimes
         write(50,'(i5,9(1x,f7.3))') time(j),Mbol(j),V(j),
     $        mag(j,43)-mag(j,44),
     $        mag(j,45)-mag(j,42),
     $        mag(j,42)-mag(j,11),
     $        mag(j,42)-mag(j,5),
     $        mag(j,42)-mag(j,6),
     $        mag(j,9)-mag(j,10),
     $        mag(j,10)-mag(j,11)
      end do
      write(50,55)
     $     'time','K-L','L-M','V-RJ','V-IJ',
     $     'JK-V','UK-JK','JK-FK','FK-NK','2000-V'
      do j=1,ntimes
         write(50,'(i5,1x,9(f7.3,1x))') time(j),
     $        mag(j,11)-mag(j,12),
     $        mag(j,12)-mag(j,13),
     $        mag(j,42)-mag(j,7),
     $        mag(j,42)-mag(j,8),
     $        mag(j,15)-mag(j,42),
     $        mag(j,14)-mag(j,15),
     $        mag(j,15)-mag(j,16),
     $        mag(j,16)-mag(j,17),
     $        mag(j,18)-mag(j,42)
      end do
      write(50,56)
     $     'time',
     $     'V-ID','ID-JD','JD-KD','BJ-V','BJ-RF',
     $     'V-606','300-450','450-606','606-814'
      do j=1,ntimes
         write(50,'(i5,1x,9(f7.3,1x))') time(j),
     $        mag(j,42)-mag(j,19),
     $        mag(j,19)-mag(j,20),
     $        mag(j,20)-mag(j,21),
     $        mag(j,22)-mag(j,42),
     $        mag(j,22)-mag(j,23),
     $        mag(j,42)-mag(j,26),
     $        mag(j,24)-mag(j,25),
     $        mag(j,25)-mag(j,26),
     $        mag(j,26)-mag(j,27)
      end do
      write(50,57)
     $     'time','u''-g''','g''-r''','V-r''',
     $     'r''-i''','i''-z''','u-v','v-g','g-V','g-r'
      do j=1,ntimes
         write(50,'(i5,1x,9(f7.3,1x))') time(j),
     $        mag(j,28)-mag(j,29),
     $        mag(j,29)-mag(j,30),
     $        mag(j,42)-mag(j,30),
     $        mag(j,30)-mag(j,31),
     $        mag(j,31)-mag(j,32),
     $        mag(j,33)-mag(j,34),
     $        mag(j,34)-mag(j,35),
     $        mag(j,35)-mag(j,42),
     $        mag(j,35)-mag(j,36)
      end do
      write(50,58) 
     $     'time','1650-B','1650-2500','3150-B'
      do j=1,ntimes
         write(50,'(i5,3(1x,f7.3))') time(j),
     $        mag(j,39)-mag(j,45),
     $        mag(j,39)-mag(j,40),
     $        mag(j,41)-mag(j,45)
      end do
      close(50)

 51   format(t2,a,t9,a,t19,a,t28,a,t36,a,t45,a,t54,a,t63,a,t70,a,t79,a)
 52   format(t2,a,t9,a,t18,a,t24,a,t37,a,t45,a,t54,a,t61,a,t70,a)
 53   format(t2,a,t7,a,t18,a,t27,a,t36,a,t45,a,t52,a,t61,a,t72,a)
 54   format(t2,a,t9,a,t19,a,t26,a,t34,a,t42,a,t50,a,t58,a,t66,a,t74,a)
 55   format(t2,a,t10,a,t18,a,t26,a,t34,a,t42,a,t49,a,t57,a,t65,a,t72,a)
 56   format(t2,a,t10,a,t17,a,t25,a,t34,a,t41,a,t49,a,t56,a,t64,a,t72,a)
 57   format(t2,a,t9,a,t17,a,t26,a,t33,a,t41,a,t50,a,t58,a,t66,a,t74,a)
 58   format(t2,a,t8,a,t15,a,t25,a)

      end


c--------------------------------------------------------------------------
***** Calculation of the flux in the ith filter

      subroutine calculflux(j,nfilters,nlambdafilter,
     $     lambdafilter,trans,nlambda,lambda,
     $     flux,nlines,lambdaline,fluxline,
     $     area,fluxfilter,typetrans)

* Returns fluxfilter in erg/s, assuming input fluxes are in solar units (/A)

      implicit none
      include 'peg_include.f'

      integer typetrans(nmaxfilters)
      integer nfilters,i,p,k,nlambdafilter(nmaxfilters),nlines,nlambda,j
      double precision fluxtrans,lambdafilter(nmaxfilters,1000)
      double precision trans(nmaxfilters,1000),lambda(nmaxlambda)
      double precision area(nmaxfilters),flux(nmaxlambda)
      double precision lambdaline(nmaxlines),fluxline(nmaxlines)
      double precision fluxfilter(nmaxotimes,nmaxfilters)
      double precision lambdainf,lambdasup,transline
      double precision fluxinf,fluxsup,transinf,transsup
      double precision interplinlin,interploglog
      double precision Lsol
      parameter(Lsol=3.826e33)

      do i=1,nfilters
         fluxtrans=0.
         if ((lambda(1).gt.lambdafilter(i,1)).or.(lambda(nlambda).lt.
     $        lambdafilter(i,nlambdafilter(i)))) then
            fluxtrans=-1.
         else
            k=1
            do while(lambda(k+1).lt.lambdafilter(i,1))
               k=k+1
            enddo
            p=1
            lambdainf=lambdafilter(i,p)
            transinf=0.
            fluxinf=0.
            lambdasup=lambdainf
            do while(p+1.le.nlambdafilter(i).and.k+1.le.nlambda)
               if (lambdafilter(i,p+1).lt.lambda(k+1)) then
                  lambdasup=lambdafilter(i,p+1)
                  fluxsup=interploglog(lambda(k),lambda(k+1),flux(k),
     $                 flux(k+1),lambdasup)
                  transsup=trans(i,p+1)*lambdasup**typetrans(i)         
                  p=p+1
               else
                  lambdasup=lambda(k+1)
                  fluxsup=flux(k+1)
                  transsup=interplinlin(lambdafilter(i,p),
     $                 lambdafilter(i,p+1),trans(i,p),trans(i,p+1),
     $                 lambdasup)*lambdasup**typetrans(i)  
                  k=k+1
               endif
               fluxtrans=fluxtrans+(transinf*fluxinf
     $              +transsup*fluxsup)*(lambdasup-lambdainf)/2
               lambdainf=lambdasup
               fluxinf=fluxsup
               transinf=transsup
            enddo
            do k=1,nlines
               p=1
               if ((lambdaline(k).gt.lambdafilter(i,1)).and.
     $              (lambdaline(k).lt.lambdafilter(i,nlambdafilter(i))))
     $              then
                  do while (.not.((lambdaline(k).ge.lambdafilter(i,p))
     $                 .and.(lambdaline(k).lt.lambdafilter(i,p+1))))
                     p=p+1
                  end do
                  transline=(trans(i,p)
     $                 +(lambdaline(k)-lambdafilter(i,p)) 
     $                 *(trans(i,p+1)-trans(i,p))
     $                 /(lambdafilter(i,p+1)-lambdafilter(i,p))) 
     $                 *lambdaline(k)**typetrans(i)
                  fluxtrans=fluxtrans+transline*fluxline(k)
               end if
            end do           
c           AL- from Lsun to erg/s:
            fluxfilter(j,i)=fluxtrans*Lsol/area(i)
         endif
      end do
 
      end
      
***** 

      double precision function interplinlin(x1,x2,y1,y2,t)

      implicit none
      double precision x1,x2,y1,y2,t

      interplinlin=y1+(y2-y1)*(t-x1)/(x2-x1)

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

      end

