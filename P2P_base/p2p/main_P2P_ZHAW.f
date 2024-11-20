      program P2P_REA

C Dr. Frank Holzaepfel, Stephan Koerner, Institute of Atmospheric Physics, German Aerospace Center,
C Germany, Oberpfaffenhofen
C 10.06.2014

C This code and the respective object files are provided to NASA in the scope of 
C a NASA-DLR cooperation which comprises among others the exchange of wake vortex models.
C Both this code and the object files may only be used according to the conditions described 
C in the software exchange agreement. 

C shear gradient interpolated, apparent gradient of shear avoided by extrapolation
C of lateral wind (nzmax), depending on sh2 renamed to dgamsh

C geometrically similar compression of distance between primary & sec. vortices
C for vortex generation at z < 0.95 hsec (minimum generation height: 0.1 b0)

C deterministic parametrization of wind shear and shear gradient for Sodar lateral winds
C developed on: 04.07.2006, implemented on: 14.12.2006

C general idea: circulation of wv decreases in 2 phases; radii averaged circulation (e.g. Gamma_3-8) 
C is to be used model is dimensionless

C Phase 1: internal and turbulent diffusion with effective viscosity n1=0.00178
C initial structure of the fully developed vortice is determined by its age T1=-2.22
C and the constant A=1.1

C Phase 2: fast decay triggered by different mechanisms of instability starts
C at time T2 with effective viscosity n2;
C T2, nu2 =f(N*,EDR*,Sh,..)

C If stratification is stable w is additionally decelerated according to Greene

C case id from casefile
C axial, vertical and lateral winds
C effects of headwind & flight path angle on vertical transport and vortex age considered in fixed prediction plane

C height correction for axial wind for z < himage = 0.6 b0 deactivated

C correction of Gamma_5-15 output (Gam0=Gam0*Gamfac)

C T2=1.2*T2 allows for long-lasting A340 vortices

C increased safety factor for sh*>1 (06.11.03)

C Cydsh, Czdsh = 0.3 (23.1.04)

C axial wind shear implemented (quadratic superposition)
C influence of shear limited (dsh=min(dsh,Gam/Czdsh))
C Cydsh = 0.3, Czdsh = 0.4 (10.2.04)



      implicit none

      integer*8 nn
      parameter (nn=250*25) ! max vortex age 250 s

      integer*8 nndit, i, j, k, kk, l, nit, it, icase, il, GE, p
      integer*8 ho, hi, nz, nEDR, nnz, nzmax, nnedr
      integer*8 ny, nc, ntp, nts
      integer*8 geinit1(2),geinit2(2),geinit3(2)
      integer*8 ifirst,ih,hour,im,mi,idum,acID,waitGE,sec
      integer*8 iGtab,nlid,cs_lidar,klid,klidl,timeshift
      integer*8 cs_um,cs_vm,empty_lines,col
      integer*8 lines_meteo,EDR_on
      integer*8 d_t(8),tms,seed(3)
      integer*8 normalize,r_lo,r_up,r_auto

      PARAMETER (nndit=20,nnz=200,nlid=200)
      PARAMETER (nnedr=200)		

      real*8 t, y(2), x0, y0, z(2), z0, zweight, dz, gpa, r
      real*8 b0, Gam(2), Gam0, Gamold(2), Gamfac, Gammin, Gammax
      real*8 N, Nmean, sh, dsh, shsq, dtsh, CshW, sh2, sh21, sh22,Csh2
      real*8 dt, t0, w0
      real*8 A, nu1, nu2(2), T1, T2
      real*8 pi, C, rho
      real*8 weight,temperature(nnz),static_pressure(nnz)
      real*8 Gamtab(120),rctab(120),G515
      real*8 Gamb
      real*8 qq(nnz),tke(nnz)
      real*8 ydel,zdelp,zdelm,Cydel,Cydel0,Czdel,Czdel0,Cydsh,Czdsh,db
      real*8 d,dum,gami
      real*8 utemp,umean,uac,dtu,dtufac,dzu,zmean
      real*8 vtemp(2),wtemp,qtemp,dtemp,thetemp
      real*8 vGE
      real*8 EDR,EDRmean,EDR90,EDR60,EDRSF
      real*8 GPSt,umag,udir,upmean,uqmean,up,uq,head,BS,dtA
      real*8 uqrms,qATTAS
      real*8 um(nnz),vm(nnz)
      real*8 zm(nnz),wm(nnz),tm(nnz),sw(nnz)
      real*8 th(nnz),thref,vtraj
      real*8 zEDR(nnedr),EDRpr(nnedr)
      real*8 out(2,6,nn),out_mean(7,nn)
      real*8 tp,yp,zp,cp,ts,ys,zs,cs,rmsy,rmsz,rmsc,tiny,tdif,tdifold
      real*8 y1,z1,y2,z2
      real*8 y3, typ, tzp, tym, tzm
      real*8 tmax
      real*8 ynd,y1nd,y2nd,ysto(13,2,2),prob(13),ypuJo,ysuJo
      real*8 znd,zsto(13,2,2),zuJo,gnd,gsto(13,2),guJo
      real*8 ylo,yup,zlo,zup,glo,gup
      real*8 yloc(12),zloc(12)
      real*8 zb0,Tg,TgP
      real*8 mass,span,gamlf,gamrf,tfirst,tlf,trf,ylf,yrf,zlf,zrf
      real*8 gamls,gamrs,gamlt,gamrt
      real*8 gaml(nlid),gamr(nlid),tl(nlid),tr(nlid),yl(nlid),yr(nlid)
      real*8 zl(nlid),zr(nlid)
      real*8 tlid, tDFS, tdiff, Gam0DLH
      real*8 t00,Gamlin(2),Gam00(2),zold
      real*8 edrson,qual,shmax,sh2max,dgamsh
      real*8 tfac,nu2fac
      real*8 z0_dis(100),wind_dir_lid
      real*8 wind_dir,b(nlid),EDR_scatter(100),EDR_z0,zm_scatter(100)
      real*8 biasyr(10000),biaszr(10000),biascr(10000),tbiasr(10000)
      real*8 biasyl(10000),biaszl(10000),biascl(10000),tbiasl(10000)
      real*8 biasb(10000)
      real*8 nan_met,nan_lid
      real*8 height_offset_lidar,lat_offset_lidar
      real*8 inboard_loading,weakening(2),thin1(2)


      namelist /DIR/ filename

      character*100 filename, mfilename, ofilename,output
      character*60  case,title,scordat                        
      character*8 time
      character*6 cdum,date
      character*4 ac,rwy,ac_type(100)
      character*100 ERGNAME
      character*15 b0_str
      character*3 meteo_structure(20)
      character*3 camp,rw
      character*35 edrfile,acfile
      character*21 case_new


      logical ACK
      LOGICAL :: bias_exists

      common /tab/rctab,Gamtab
      common /geinit/geinit1,geinit2,geinit3,yloc,zloc

      parameter (normalize=0) 		! if normalize=1 the output will be normalized, if 0 output will be dimensional

      character(len=128) :: ac_filename, meteo_filename, edr_filename
      character(len=500) :: input_dir, res_dir, arg, abs_path, path_name
      character(len=20) :: runpath
      character(len=100) :: run_id
      integer :: i_iteration, argcount
      logical :: abs_path_given, verbose
      ! Set default values
      input_dir = 'inputs/'
      res_dir = 'results/'
      abs_path_given=.FALSE.
      verbose=.FALSE.

      ! Process command line arguments if provided
      argcount = command_argument_count()
      if (argcount > 1) then
            ! do i_iteration = 1, argcount
            i_iteration = 1
            do while (i_iteration <= argcount)
                  call get_command_argument(i_iteration, arg)
                  if (arg == '-o') then
                        call get_command_argument(i_iteration+1, arg)
                        read(arg, *) run_id
                        i_iteration = i_iteration+1
                  endif
                  if (arg == '-p') then
                        call get_command_argument(i_iteration+1, arg)
                        abs_path = arg 
                        i_iteration = i_iteration+1
                        abs_path_given=.TRUE.
                  endif
                  if (arg == '-v') then
                        verbose=.TRUE.
                  endif
                  i_iteration = i_iteration+1
            end do
      endif

      input_dir = trim(input_dir)
      res_dir = trim(res_dir)
      abs_path = trim(abs_path)
      if (abs_path_given) then
        CALL PrintIfVerbose('Abs path: ' // trim(abs_path), verbose)
        input_dir = trim(abs_path) // '/' // input_dir // '/'
        res_dir = trim(abs_path) // '/' // res_dir // '/'
      endif

      if (len_trim(run_id) > 0) then
            CALL PrintIfVerbose('Run ID: ' // trim(run_id), verbose)
            runpath =  '/' // trim(run_id) // '/'
            ! input_dir = input_dir // '/' // trim(run_id) // '/' // res_dir
            input_dir = trim(input_dir) // trim(run_id)
            res_dir = trim(res_dir) // trim(run_id)
      endif

      CALL PrintIfVerbose('Input dir: ' // trim(input_dir), verbose)
      CALL PrintIfVerbose('Results dir: ' // trim(res_dir), verbose)

      pi=4.*atan(1.)
      nit=350      !previously 250
      gamfac=1./0.9582
      tiny=0.000001
      zb0=1.
      t00=999.
      ACK=.TRUE.    ! true if aircraft type is unknown
      EDR_on=0
      nan_met=0


c     model constants
      C=0.0121
      A=1.1
      T1=-3.48       ! Gam(0)=0.9582
c      T1=-3.        ! Gam(0)=1 (exakt 0.99626)
      nu1=0.00178
      db=1.0
      Cydel=1.0
      Czdel=0.5
      Czdel0=Czdel
      Cydel0=Cydel
      Cydsh=0.3
      Czdsh=0.4
c     for implementation of det. shear parameterization replace z(1) resp z(2) by (z(1)+z(2))/2.
      CshW=0.05    ! for sodar shear / Reduzierung der Sinkgeschw. in Scherschichten f�r beide Wirbel
      Csh2=0.2   ! for sodar shear / obiges f�r Einzelw. als Fkt. von 2ter Ableitung des CW-profils
c       CshW=0.2    ! for lidar shear / Reduzierung der Sinkgeschw. in Scherschichten f�r beide Wirbel
c       Csh2=0.3     ! for lidar shear / obiges f�r Einzelw. als Fkt. von 2ter Ableitung des CW-profils
      nzmax=nnz

c     read relevant case data used for initialization, date, time from file fort.13
      path_name = trim(input_dir) // '/fort.13'
      CALL PrintIfVerbose('Reading file ' // trim(path_name), verbose)
      open(unit=13, file=path_name)   
      read (13,*) case
      read (13,*) ifirst
      read (13,*) ac
      read (13,*) date
      read (13,*) hour
      read (13,*) mi
      read (13,*) sec
      read (13,*) camp
      read (13,*) rw
      read (13,*) edrfile
      read (13,*) acfile
      close(13)

      if(ifirst.eq.0) goto 993

c     read aircraft initial conditions from ac_init.dat
      path_name = trim(input_dir) // '/ac_init.dat'
      CALL PrintIfVerbose('Reading file ' // trim(path_name), verbose)
      open(14, file=path_name, status='old', err=998)
      read (14,*)
      read (14,*,end=19) z0,gam0,b0,uac,x0,y0,gpa	! gpa = glide path angle
19    continue
      close (14)

      !     gam0 aus LFZ Parametern (mass: aircraft mass in kg, rho: air density; b0: initial vortex separation)

!     b0 = pi/4 * span
!     gam0=(9.81*mass)/(rho*b0*uac)  

      gpa=gpa*0.01745329    ! deg to rad


c     read lidar data
      
      open (15,file='lidar.dat',status='old',err=998)
      klid=1
      do k = 1,nlid
       read (15,*,end=30)
     &      tr(k),zr(k),yr(k),gamr(k),tl(k),zl(k),yl(k),gaml(k)

       klid=klid+1
       if(k.eq.ifirst) then
!        tfirst=tl(k)
	ylf=yl(k)
	zlf=zl(k)
	gamlf=gaml(k)
	yrf=yr(k)
	zrf=zr(k)
	gamrf=gamr(k)
	trf=tr(k)
	tlf=tl(k)
        tfirst=(tlf+trf)/2 !
       endif
       if(k.eq.2) then
	gamls=gaml(k)
	gamrs=gamr(k)
       endif
       if(k.eq.3) then
	gamlt=gaml(k)
	gamrt=gamr(k)
       endif
      enddo
30    continue
      klidl=k-1
      close (15)
      
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      t0=2.*pi*b0**2/gam0              !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!     !!!!!!!!!!!!!!!!!!!!!
      tfirst=0.           !           time for which P2P starts calculating values
!      tfirst=(tlf+trf)/2 !
!     !!!!!!!!!!!!!!!!!!!!!

      if(tfirst.le.0.) tfirst=0.

      tfirst=tfirst/t0

      open (15,file='lidar_norm.dat')
      do k = 1,klidl
       write(15,990)tr(k)/t0,zr(k)/b0,yr(k)/b0,
     &      gamr(k)/gam0*(A-exp(-C/(nu1*(tfirst-T1))))*gamfac,
     &      tl(k)/t0,zl(k)/b0,yl(k)/b0,
     &      gaml(k)/gam0*(A-exp(-C/(nu1*(tfirst-T1))))*gamfac
      enddo
      close (15)

     
      write(*,*)
      write(*,'(A7,F7.1)') 't0    =   ',t0
      write(*,'(A7,F7.1)') 'gam0  =   ',Gam0
      write(*,'(A7,F7.1)')  'b0    =  ',b0
      write(*,'(A7,F7.1)')  'y0    =  ',y0
      write(*,'(A7,F7.1)')  'x0    =  ',x0
      write(*,'(A7,F7.1)')  'z0    =  ',z0
      write(*,'(A7,I7)')    'nit   =  ',nit
      write(*,'(A7,F7.1)')  'uac   =  ',uac
      write(*,'(A8,F6.1)')  'gpa   = ',gpa/0.01745329
 

c     rms deviation of predictions and measurements


      path_name = trim(res_dir) // '/RMS.dat'
      OPEN(16,FILE=path_name,access='append')


C     probability density distribution

!       OPEN(24,FILE='pddy.dat',access='append')
!       OPEN(25,FILE='pddz.dat',access='append')
!       OPEN(26,FILE='pddG.dat',access='append')



!c     create table for descent speeds for w=f(rc(Gamma_5-15)) or w=f(rc(Gamma_5-15))
!c     (original P2P version described in JoA 2003 for medium and heavy a/c) 
!
!
!c      create table for descent speeds for w=f(rc(Gamma_5-15))  
!
!       if (b0.gt.30.0) then
!        print *,'gamma average radii: 5m-15m'
!        do j=1,120
!          rctab(j) = 121.0 - j
!          G515 = 0.0
!          do k=5,15
!            G515 = G515 + 1. - exp(-1.257 * k**2 / rctab(j)**2)
!          enddo
!          Gamtab(j)=G515/11.0
!        enddo
!
!       else
!
!c      create table for descent speeds for w=f(rc(Gamma_3-8))  
!
!        print *,'gamma average radii: 3m-8m'
!        do j=1,120
!          rctab(j) = 121.0 - j
!          G515 = 0.0
!          do k=3,8
!            G515 = G515 + 1. - exp(-1.257 * k**2 / rctab(j)**2)
!          enddo
!          Gamtab(j)=G515/6.
!       endif
!      endif

c     descent speeds for consistently nondimensional model
c     create table for descent speeds for w=f(rc(Gamma_0.21*b0/2_0.65*b0/2)), corresponds to Gam5-15 for a/c with span=60m

      do j=1,120
        rctab(j) = 121.0 - j
        G515 = 0.0
        do k=0,10
            r=(0.21+k*0.044)*b0/2
            G515 = G515 + 1. - exp(-1.257 * r**2 / rctab(j)**2)
        enddo
        Gamtab(j)=G515/11.
      enddo

C     result files
      path_name = trim(res_dir) // '/TRAJEC.dat'
      OPEN(10,FILE=path_name,STATUS='UNKNOWN')

!       ERGNAME='TRAJEClp.DAT'   ! lower envelope
! 
!       OPEN(11,FILE=ERGNAME,STATUS='UNKNOWN')
! 
!       ERGNAME='TRAJECus.DAT'   ! upper envelope
! 
!       OPEN(12,FILE=ERGNAME,STATUS='UNKNOWN')

      path_name = trim(res_dir) // '/P2P_lev_y.dat'
      OPEN(20,FILE=path_name,STATUS='UNKNOWN')

      path_name = trim(res_dir) // '/P2P_lev_z.dat'
      OPEN(21,FILE=path_name,STATUS='UNKNOWN')

      path_name = trim(res_dir) // '/P2P_lev_g.dat'
      OPEN(22,FILE=path_name,STATUS='UNKNOWN')

!       ERGNAME='P2P_0-sigma.dat' ! 50th percentile 
! 
!       OPEN(23,FILE=ERGNAME,STATUS='UNKNOWN')


c     derived quantities and normalization

      w0=Gam0/(2.*pi*b0)
      y0=y0/b0
      dt=1./t0
      uac=uac/w0
      gami=Gam0

      do i=1,nn
       do j=1,6
        out(1,j,i)=9999.             ! initialize output variables
        out(2,j,i)=-9999. 
       enddo
      enddo


c     read the meteorological profiles
      path_name = trim(input_dir) // '/meteo.dat'
      open (14,file=path_name,status='old',err=998)
      read(14,*)
      do k =1,nnz
        read (14,*,end=9) zm(k),um(k),vm(k),wm(k),qq(k),
     &                    temperature(k),th(k),tke(k)
      enddo
9     continue
      close(14)


c     count number of lines in meteo.dat

      lines_meteo=0
      open (14,file=path_name,status='old',err=998)
      read(14,*)
      do k=1,nnz
         read (14,*,end=11)  
         lines_meteo=lines_meteo+1
      enddo	    
11    continue
      close(14)
      
      
c     Calculation of EDR from TKE according to Donaldson & Bilanin


c      if (edrfile.eq.'none') then
c       do k = 1,nnz
c         if (tke(k).ne.nan_met) then                  
c           qq(k)=(2.*tke(k))**0.5			
cc          qq(k)=(2./3*tke(k))**0.5			
c         endif
c       enddo
c      endif


c      if (edrfile.eq.'none') then
c      do k = 1,nnz
c         if(zm(k).lt.169.) then
c           EDRpr(k)=qq(k)**3/(8.*0.65*zm(k))
c         else
c           EDRpr(k)=qq(k)**3/(8.*110)
c         endif
c      enddo
c      endif


c     Read EDR-file if it exists

c      if (edrfile.ne.'none') then
      path_name = trim(input_dir) // '/edr.dat'
      open (14,file= path_name,status='old',err=998)
        do k=1,nnz
         read (14,*,end=21) dum,EDRpr(k) 
        enddo 
21      continue
        close (14)
c      endif

c     determine meteo index at initial height of the vortex pair     

      do kk=1,nnz
         if (zm(kk).ge.z0) then
c         nz=kk+10                   ! meteo 10 levels above initial height to cover rebound
          nz=kk+5                    ! meteo  5 levels above initial height to cover rebound
c         nz=lines_meteo		     ! meteo according to number of measured heights in meteo.dat
         goto 24
         endif
      enddo
      write(*,*)'meteo data missing'
      goto 999
24    continue


c     determine EDR index at initial height of the vortex pair
      nedr=nnz
      do kk=1,nedr-1
c         write(*,*) kk,zEDR(kk),z0
         if (zEDR(kk).ge.z0) then
         nEDR=kk+4                  ! meteo 4 levels above initial height to cover rebound
          goto 12
         endif
      enddo

12    continue


C     normalize z0

      z0=z0/b0


C     parameters for z0 < b0

      z0=max(0.1,z0)
      tfac=z0
      if(z0.lt.0.125) tfac=tfac*(z0/0.125)**2
      tfac=min(1.,tfac)
      nit=int(nit/tfac)
      dt=1./t0*tfac

      nu2fac=1.
      if(z0.le.1.) then
        nu2fac=0.4/(0.00473+1.13*z0-1.05*z0**2+0.323*z0**3)
        nu2fac=nu2fac**2
      endif



c     inter- and extrapolation of met data

      do k = 1,nz+1
       if(vm(k).eq.nan_met) then
          nzmax=k-1
	  goto 1012
       endif
      enddo
1012  continue


      do k = 2,nz+1
       if(abs(um(1)).eq.nan_met) um(1)=um(2)
       if(abs(um(k)).eq.nan_met) then
         um(k)=(um(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +um(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(abs(um(k-1)).eq.nan_met.or.abs(um(k+1)).eq.nan_met) then
          um(k)=um(k-1)
         endif
       endif
       if(abs(vm(1)).eq.nan_met) vm(1)=vm(2)
       if(abs(vm(k)).eq.nan_met) then
         vm(k)=(vm(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +vm(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(abs(vm(k-1)).eq.nan_met.or.abs(vm(k+1)).eq.nan_met) then
          vm(k)=vm(k-1)
         endif
       endif
      enddo

      do k = 2,nz
       if(sw(1).eq.nan_met) sw(1)=sw(2)
       if(sw(k).eq.nan_met) then
         sw(k)=(sw(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +sw(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(sw(k-1).eq.nan_met.or.sw(k+1).eq.nan_met) then
          sw(k)=sw(k-1)
         endif
       endif
       if(wm(1).eq.nan_met) wm(1)=wm(2)
       if(wm(k).eq.nan_met) then
         wm(k)=(wm(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +wm(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(wm(k-1).eq.nan_met.or.wm(k+1).eq.nan_met) then
          wm(k)=wm(k-1)
         endif
       endif
       if(qq(1).eq.nan_met) qq(1)=qq(2)
       if(qq(k).eq.nan_met) then
         qq(k)=(qq(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +qq(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(qq(k-1).eq.nan_met.or.qq(k+1).eq.nan_met) then
          qq(k)=qq(k-1)
         endif
       endif
       if(th(1).eq.nan_met) th(1)=th(2)
       if(th(k).eq.nan_met) then
         th(k)=(th(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +th(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(th(k+1).eq.nan_met) then
c      Extrapolation of temperature (N=0)
         th(k)=(zm(k)-zm(k-1))/(zm(k-1)-zm(k-2))*(th(k-1)
     &         -th(k-2))+th(k-1)
c         th(k)=th(k-1)
c      Extrapolation of temperature gradient
c          th(k)=th(k-1)
c     &    +(th(k-1)-th(k-2))*(zm(k)-zm(k-1))/(zm(k-1)-zm(k-2))
c      halb/halb
c           th(k)=th(k-1)
c           th(k)=0.5*(th(k)+th(k-1)
c      &    +(th(k-1)-th(k-2))*(zm(k)-zm(k-1))/(zm(k-1)-zm(k-2)))
         endif
       endif
       if((EDRpr(1).eq.nan_met).or.(EDRpr(1).eq.0)) EDRpr(1)=EDRpr(2)
       if(EDRpr(k).eq.nan_met) then
         EDRpr(k)=(EDRpr(k-1)*(zm(k+1)-zm(k))/(zm(k+1)-zm(k-1))
     &         +EDRpr(k+1)*(zm(k)-zm(k-1))/(zm(k+1)-zm(k-1)))
         if(EDRpr(k-1).eq.nan_met.or.EDRpr(k+1).eq.nan_met) then
          EDRpr(k)=EDRpr(k-1)
         endif
       endif
      enddo


c     normalize the meteorological data:

      do i=1,nz
         zm(i)=zm(i)/b0
         qq(i)=qq(i)/w0
         um(i)=um(i)/w0
         vm(i)=vm(i)/w0
         wm(i)=wm(i)/w0
         EDRpr(i)=(EDRpr(i)*b0)**(1./3.)/w0
      enddo


c     determine the normalized initial height of the vortex pair

      do kk=2,nnz
         if (zm(kk).ge.z0) then
         ho=kk-1
         goto 25
         endif
      enddo
      goto 999
25    continue

      zweight=(z0-zm(ho))/(zm(ho+1)-zm(ho))
      thref=(1.-zweight)*th(ho)+zweight*th(ho+1)


C     crosswind direction (according to measurements)

      if (((yr(4)+yl(4))/2.lt.(yr(1)+yl(1))/2).and. 
     &    ((yr(3)+yl(3))/2.lt.(yr(1)+yl(1))/2).and.
     &    ((yr(2)+yl(2))/2.lt.(yr(1)+yl(1))/2)) then        	! use of vortex drift instead of wind data for post-processing 
         wind_dir_lid=-1			        	
      else
         wind_dir_lid=1
      endif



c     3 runs with different T2 and nu2 

      do 50 icase=0,2



c     initialization

      GE=0
      y(1)=y0-0.5
      y(2)=y0+0.5
      z=z0
      ydel=0.
      zdelp=0.
      zdelm=0.
      dsh=0.
      dtsh=0.
      dgamsh=0.
      Gam=1.
      Gamb=1.             
      T2=999.              
      Tg=999.
      tgP=999.
      waitGE=0
      EDRmean = 0.
      Nmean = 0.
      umean = 0.
      iGtab=120
      weakening=1
      thin1=1.75001*pi  ! rotation angle angle of sec. vortices. Initial angle = 315 �

c     this is the big loop

      do 100 it=1,nit
      t=tfirst+it*dt


c     determine the height of the vortex pair in order to determine the meteorological parameters

      do kk=2,nnz
        if (zm(kk).ge.z(1)) then
         hi=kk-1
         goto 261
        endif
      enddo

261   continue

      zweight=(z(1)-zm(hi))/(zm(hi+1)-zm(hi))
      vtemp(1)=(1.-zweight)*vm(hi)+zweight*vm(hi+1)

       do kk=2,nnz
         if (zm(kk).ge.z(2)) then
          hi=kk-1
          goto 262
         endif
       enddo

262    continue

      zweight=(z(2)-zm(hi))/(zm(hi+1)-zm(hi))
      vtemp(2)=(1.-zweight)*vm(hi)+zweight*vm(hi+1)

       zmean=(z(1)+z(2))/2.

       do kk=2,nnz
         if (zm(kk).ge.zmean) then
          hi=kk-1
          goto 26
         endif
       enddo

26    continue

      zweight=(zmean-zm(hi))/(zm(hi+1)-zm(hi))
      utemp=(1.-zweight)*um(hi)+zweight*um(hi+1)

c compute crosswind shear, total shear, and derivative of shear of Sodar data

      sh = (((vm(hi+1)-vm(hi))/(zm(hi+1)-zm(hi)))**2
     &    +  ((um(hi+1)-um(hi))/(zm(hi+1)-zm(hi)))**2)**0.5
c       sh  = (vm(hi+1)-vm(hi))/(zm(hi+1)-zm(hi))

	if (hi.eq.1) then
	  sh2 = 0.
	elseif (hi+2.gt.nzmax) then
	  sh2 = 0.
	else
      	  sh22 = (vm(hi+2)-2.*vm(hi+1)+vm(hi))/
     &           (0.5*(zm(hi+2)-zm(hi)))**2
      	  sh21 = (vm(hi+1)-2.*vm(hi)+vm(hi-1))/
     &           (0.5*(zm(hi+1)-zm(hi-1)))**2
          zweight = (zmean-zm(hi))/(zm(hi+1)-zm(hi))
          sh2 = (1.-zweight)*sh21+zweight*sh22
    	endif
	sh2=0.

       dsh=dsh*(1.-exp(-5./(t-dtsh)))
       if(sh.gt.1.0.and.sh.gt.dsh) then
         dsh=sh
         dtsh=t
       endif

       shsq=((um(hi+1)-um(hi))/(zm(hi+1)-zm(hi)))**2
     &     +((vm(hi+1)-vm(hi))/(zm(hi+1)-zm(hi)))**2


c determine the height of the vortex pair in order to determine the meteorological parameters

       do kk=2,nnz
         if (zm(kk).ge.zmean) then
          hi=kk-1
          goto 28
         endif
       enddo

28    continue

      zweight=(zmean-zm(hi))/(zm(hi+1)-zm(hi))
      qtemp=(1.-zweight)*qq(hi)+zweight*qq(hi+1)
      wtemp=(1.-zweight)*wm(hi)+zweight*wm(hi+1)
      thetemp=(1.-zweight)*th(hi)+zweight*th(hi+1)


c     compute the Brunt-Vaisala frequency N for WV
c     relative to the initial height of the WV

      dtemp=thref-thetemp
      N=0.
      if(dtemp/(z0-(zmean)+tiny).gt.0.0) then
        N=(9.81/thref*dtemp/((z0-(zmean)+tiny)*b0))**0.5*t0
      endif


!      EDR data
       EDR=(1.-zweight)*EDRpr(hi)+zweight*EDRpr(hi+1)

c     utemp=0.
      wtemp=0.

      EDRmean=((it-1)*EDRmean+EDR)/it  ! running average
      Nmean=((it-1)*Nmean+N)/it
      umean=((it-1)*umean+utemp)/it

c      print *,Nmean,'Nmean'
c      print *,EDRmean,'EDRmean'

      dtu=(umean/(uac+umean))/(1-umean/(uac+umean))  ! ground speed = uac+umean, dtu=timeshift due to headwind
      dzu=tan(gpa)*(uac+umean)*dt*dtu
      dtu=0.
      dzu=0.

c      activation of in-ground-decay: Tg=t(z=1)

       if(t.ge.T2) goto 14
       if(t.ge.Tg) goto 14

       if(zmean.le.zb0*1.0) then

	Tg=t

c      determine crosswind at z=0.6 b0 for IGE decay

       do kk=2,nnz
         if (zm(kk).ge.0.6) then
          hi=kk-1
          goto 290
         endif
       enddo
290    continue

      zweight=(0.6-zm(hi))/(zm(hi+1)-zm(hi))
      vGE=(1.-zweight)*vm(hi)+zweight*vm(hi+1)

      endif

14    continue

       Gamold=Gam

       Gam=A-exp(-C/(nu1*(t-T1)))

       call T2_nu2 (icase,Nmean,EDRmean,EDRpr(1),vGE,Tg,T2,
     &              nu2,nu2fac,waitGE)

       if (t.ge.T2) then
        Gam=Gam-exp(-C/(nu2*(t-T2)))
       endif

       if (Gam(1).le.-0.6) goto 101


       call traj (Gam,Gamb,b0,y,z0,z,dt,N,vtemp,wtemp,vGE,GE,dzu,
     &      dsh,CshW,sh2,Csh2,dgamsh,iGtab,icase,weakening,thin1)


       dsh=min(dsh,abs(Gam(1)/Czdsh)) ! limitation of shear effects to descent speed

       if(zmean*b0.lt.100.) then
        Czdel=Czdel0*zmean*b0/100. ! linear damping of vertical turbulent spreading in ground proximity
        Cydel=Cydel0*(1.5-0.5*zmean*b0/100.) ! strengthen lateral turbulent spreading in ground proximity
      endif

       ydel=ydel+((qtemp*Cydel)**2+(dsh*Cydsh)**2)**0.5*dt
       if(t.le.Tg)zdelp=zdelp+((qtemp*Czdel)**2+(dsh*Czdsh)**2)**0.5*dt
       zdelm=zdelm+qtemp*dt*Czdel


c     Normalized Quantities for Output

       if(icase.eq.0) then
         out_mean(1,it)=t*(1.-dtu)                        ! out_mean describes deterministic behavior
         out_mean(2,it)=-Gam(1)*Gamfac
         out_mean(3,it)=y(1)
         out_mean(4,it)=z(1)
         out_mean(5,it)=Gam(2)*Gamfac
         out_mean(6,it)=y(2)
         out_mean(7,it)=z(2)
       endif


c     Normalized Quantities for Output

        out(1,1,it)=t*(1.-dtu)
        Gammin=min(Gam(1),Gam(2))
        if ((Gammin-0.2)*Gamfac.lt.out(1,2,it))
     &     out(1,2,it)=(Gammin-0.2)*Gamfac     
        if ((y(1)-ydel-db).lt.out(1,3,it))
     &     out(1,3,it)=(y(1)-ydel-db)
        if ((y(2)-ydel-db).lt.out(1,4,it))
     &     out(1,4,it)=(y(2)-ydel-db)
        if ((z(1)-zdelm-db).lt.out(1,5,it))
     &     out(1,5,it)=max(0.,(z(1)-zdelm-db))
        if ((z(2)-zdelm-db).lt.out(1,6,it))
     &     out(1,6,it)=max(0.,(z(2)-zdelm-db))


        out(2,1,it)=t*(1.-dtu)
        Gammax=max(Gam(1),Gam(2))
        if ((Gammax+0.2)*Gamfac.gt.out(2,2,it))
     &     out(2,2,it)=(Gammax+0.2)*Gamfac
        if ((y(1)+ydel+db).gt.out(2,3,it))
     &     out(2,3,it)=(y(1)+ydel+db)
        if ((y(2)+ydel+db).gt.out(2,4,it))
     &     out(2,4,it)=(y(2)+ydel+db)
        if ((z(1)+zdelp+db).gt.out(2,5,it))
     &     out(2,5,it)=(z(1)+zdelp+db)
        if ((z(2)+zdelp+db).gt.out(2,6,it))
     &     out(2,6,it)=(z(2)+zdelp+db)


100   enddo

101   continue


c     fill up rest - normalized Quantities for Output

      do i=it,nit
         out(1,1,i)=i*(1.-dtu)/t0
         Gammin=min(Gamold(1),Gamold(2))
         if ((Gammin-0.2)*Gamfac.lt.out(1,2,i))
     &      out(1,2,i)=max(0.001,(Gammin-0.2)*Gamfac)
         if ((y(1)-ydel-db).lt.out(1,3,i))
     &      out(1,3,i)=(y(1)-ydel-db)
          if ((y(2)-ydel-db).lt.out(1,4,i))
     &      out(1,4,i)=(y(2)-ydel-db)
         if ((z(1)-zdelm-db).lt.out(1,5,i))
     &      out(1,5,i)=max(0.,(z(1)-zdelm-db))
         if ((z(2)-zdelm-db).lt.out(1,6,i))
     &      out(1,6,i)=max(0.,(z(2)-zdelm-db))
         out(2,1,i)=i*(1.-dtu)/t0
         Gammax=max(Gamold(1),Gamold(2))
         if ((Gammax+0.2)*Gamfac.gt.out(2,2,i))
     &      out(2,2,i)=(Gammax+0.2)*Gamfac
         if ((y(1)+ydel+db).gt.out(2,3,i))
     &      out(2,3,i)=(y(1)+ydel+db)
          if ((y(2)+ydel+db).gt.out(2,4,i))
     &      out(2,4,i)=(y(2)+ydel+db)
         if ((z(1)+zdelp+db).gt.out(2,5,i))
     &      out(2,5,i)=(z(1)+zdelp+db)
          if ((z(2)+zdelp+db).gt.out(2,6,i))
     &      out(2,6,i)=(z(2)+zdelp+db)


      enddo


50    enddo

C     wind direction according to model prediction

      if ((vm(1)+vm(2)+vm(3))/3.lt.0)then
         wind_dir=-1			        	
      else
         wind_dir=1
      endif


c     mean output

c     fill array with last value


      do i=1,nit
        if(out(2,2,i).le.0.)    goto 200
        if(out_mean(5,i).le.0.) goto 199
        if (normalize.eq.1) then
          write(10,61) (out_mean(j,i),j=1,7), wind_dir                 ! write TRAJEC.DAT, true wind_dir added (though calculated again in trajec_sort)
        else	
          write(10,61) out_mean(1,i)*t0,
     &       out_mean(2,i)*Gam0,
     &       out_mean(3,i)*b0,
     &       out_mean(4,i)*b0,
     &       out_mean(5,i)*Gam0,
     &       out_mean(6,i)*b0,
     &       out_mean(7,i)*b0,
     &       wind_dir
        endif
199     continue


! c       output of uncertainty range
! 
!         if (normalize.eq.1) then
!           write(11,61) (out(1,j,i),j=1,6)
!           write(12,61) (out(2,j,i),j=1,6)
!         else
!           write(11,61) out(1,1,i)*t0,
!      &	               out(1,2,i)*Gam0,
!      &	               out(1,3,i)*b0,
!      &	               out(1,4,i)*b0,
!      &                 out(1,5,i)*b0,
!      &                 out(1,6,i)*b0
!           write(12,61) out(2,1,i)*t0,
!      &	               out(2,2,i)*Gam0,
!      &	               out(2,3,i)*b0,
!      &	               out(2,4,i)*b0,
!      &                 out(2,5,i)*b0,
!      &                 out(2,6,i)*b0
! 	endif
              

c     probability levels OGE

       ylo=out(1,3,i)
       yup=out(2,4,i)

       ysto(1,1,1)=-0.2022*(yup-ylo)+ylo     !  -4 sigma, 99.994%
c       ysto(2,1,1)=-0.04799*(yup-ylo)+ylo   !  -3 sigma, 99.73%
       ysto(2,1,1)=-0.02850*(yup-ylo)+ylo   !  -3 sigma, 99.73% sym
       ysto(3,1,1)= 0.01853*(yup-ylo)+ylo   !                99.0%
c       ysto(4,1,1)= 0.1121*(yup-ylo)+ylo     !  -2 sigma, 95.4%
       ysto(4,1,1)= 0.1481*(yup-ylo)+ylo     !  -2 sigma, 95.4% sym
       ysto(5,1,1)= 0.1704*(yup-ylo)+ylo     !                90.0%
       ysto(6,1,1)= 0.2784*(yup-ylo)+ylo     !  -1 sigma, 68.3%
       ysto(7,1,1)= 0.4509*(yup-ylo)+ylo     !   0 sigma,  0.0%
       ysto(8,1,1)= 0.6300*(yup-ylo)+ylo     !   1 sigma, 68.3%
       ysto(9,1,1)= 0.7492*(yup-ylo)+ylo     !                90.0%
c       ysto(10,1,1)= 0.8160*(yup-ylo)+ylo   !   2 sigma, 95.4%
       ysto(10,1,1)= 0.8512*(yup-ylo)+ylo   !   2 sigma, 95.4% sym
       ysto(11,1,1)=0.9271*(yup-ylo)+ylo    !                99.0%
c       ysto(12,1,1)=1.009*(yup-ylo)+ylo      !   3 sigma, 99.73%
       ysto(12,1,1)=1.0285*(yup-ylo)+ylo      !   3 sigma, 99.73% sym
       ysto(13,1,1)=1.209*(yup-ylo)+ylo      !   4 sigma, 99.994%

       zlo=min(out(1,5,i),out(1,6,i))
       zup=max(out(2,5,i),out(2,6,i))

       zsto(1,1,1)=max(0.,-0.06480*(zup-zlo)+zlo)   !  -4 sigma, 99.994%
       zsto(2,1,1)= max(0.,0.09534*(zup-zlo)+zlo)   !  -3 sigma, 99.73%
       zsto(3,1,1)= max(0.,0.14880*(zup-zlo)+zlo)   !                99.0%
       zsto(4,1,1)= max(0.,0.2150*(zup-zlo)+zlo)     !  -2 sigma, 95.4%
       zsto(5,1,1)= max(0.,0.2530*(zup-zlo)+zlo)     !                90.0%
       zsto(6,1,1)= max(0.,0.3215*(zup-zlo)+zlo)     !  -1 sigma, 68.3%
       zsto(7,1,1)= max(0.,0.4390*(zup-zlo)+zlo)     !   0 sigma,  0.0%
       zsto(8,1,1)= max(0.,0.5943*(zup-zlo)+zlo)     !   1 sigma, 68.3%
       zsto(9,1,1)= max(0.,0.7304*(zup-zlo)+zlo)     !                90.0%
       zsto(10,1,1)=max(0.,0.8116*(zup-zlo)+zlo)    !   2 sigma, 95.4%
       zsto(11,1,1)=max(0.,1.008*(zup-zlo)+zlo)      !                99.0%
       zsto(12,1,1)=max(0.,1.176*(zup-zlo)+zlo)      !   3 sigma, 99.73%
       zsto(13,1,1)=max(0.,1.735*(zup-zlo)+zlo)      !   4 sigma, 99.994%

       glo=out(1,2,i)
       gup=out(2,2,i)

       gsto(1,1)=max(0.,-0.8650*(gup-glo)+glo)    ! -4 sigma,  99.994%
       gsto(2,1)=max(0.,-0.1591*(gup-glo)+glo)   ! -3 sigma,  99.73%
       gsto(3,1)= max(0.,0.02454*(gup-glo)+glo)   !                99.0%
       gsto(4,1)= max(0.,0.2118*(gup-glo)+glo)     !  -2 sigma, 95.4%
       gsto(5,1)= max(0.,0.2975*(gup-glo)+glo)     !                90.0%
       gsto(6,1)= max(0.,0.4210*(gup-glo)+glo)     !  -1 sigma, 68.3%
       gsto(7,1)= max(0.,0.5718*(gup-glo)+glo)     !   0 sigma,  0.0%
       gsto(8,1)= max(0.,0.5943*(gup-glo)+glo)     !   1 sigma, 68.3%
       gsto(9,1)= max(0.,0.7351*(gup-glo)+glo)     !                90.0%
       gsto(10,1)=max(0.,0.8823*(gup-glo)+glo)    !   2 sigma, 95.4%
       gsto(11,1)=max(0.,0.9889*(gup-glo)+glo)    !                99.0%
       gsto(12,1)=max(0.,1.454*(gup-glo)+glo)      !   3 sigma, 99.73%
       gsto(13,1)=max(0.,2.353*(gup-glo)+glo)      !   4 sigma, 99.994%

       do j=1,13
         ysto(j,1,2)=ysto(j,1,1)    ! identical cvalues for port and starb OGE
         zsto(j,1,2)=zsto(j,1,1)
       enddo

c     probability levels IGE

       ylo=out(1,3,i)
       yup=out(2,3,i)

       ysto(1,2,1)=-1.9384*(yup-ylo)+ylo     !  -4 sigma, 99.994%
c       ysto(2,2,1)=-0.5618*(yup-ylo)+ylo   !  -3 sigma, 99.73%
       ysto(2,2,1)=-0.2794*(yup-ylo)+ylo   !  -3 sigma, 99.73% sym
       ysto(3,2,1)=-0.2461*(yup-ylo)+ylo   !   -2.58        99.0%
c       ysto(4,2,1)= 0.04577*(yup-ylo)+ylo     !  -2 sigma, 95.4% 
       ysto(4,2,1)= 0.16829*(yup-ylo)+ylo     !  -2 sigma, 95.4% sym
       ysto(5,2,1)= 0.1681*(yup-ylo)+ylo     !  -1.645        90.0%
       ysto(6,2,1)= 0.3207*(yup-ylo)+ylo     !  -1 sigma, 68.3%
       ysto(7,2,1)= 0.4604*(yup-ylo)+ylo     !   0 sigma,  0.0%
       ysto(8,2,1)= 0.5649*(yup-ylo)+ylo     !   1 sigma, 68.3%
       ysto(9,2,1)= 0.6482*(yup-ylo)+ylo     !                90.0%
c       ysto(10,2,1)= 0.7092*(yup-ylo)+ylo   !   2 sigma, 95.4% 
       ysto(10,2,1)= 0.8317*(yup-ylo)+ylo   !   2 sigma, 95.4% sym
       ysto(11,2,1)=0.8489*(yup-ylo)+ylo    !                99.0%
c       ysto(12,2,1)=0.9969*(yup-ylo)+ylo      !   3 sigma, 99.73%
       ysto(12,2,1)=1.2793*(yup-ylo)+ylo      !   3 sigma, 99.73% sym
       ysto(13,2,1)=1.6343*(yup-ylo)+ylo      !   4 sigma, 99.994%

       ylo=out(1,4,i)
       yup=out(2,4,i)

       ysto(1,2,2)=-1.9384*(yup-ylo)+ylo     !  -4 sigma, 99.994%
c       ysto(2,2,2)=-0.5618*(yup-ylo)+ylo   !  -3 sigma, 99.73%
       ysto(2,2,2)=-0.2794*(yup-ylo)+ylo   !  -3 sigma, 99.73% sym
       ysto(3,2,2)=-0.2461*(yup-ylo)+ylo   !   -2.58        99.0%
c       ysto(4,2,2)= 0.04577*(yup-ylo)+ylo     !  -2 sigma, 95.4%
       ysto(4,2,2)= 0.16829*(yup-ylo)+ylo     !  -2 sigma, 95.4% sym
       ysto(5,2,2)= 0.1681*(yup-ylo)+ylo     !  -1.645        90.0%
       ysto(6,2,2)= 0.3207*(yup-ylo)+ylo     !  -1 sigma, 68.3%
       ysto(7,2,2)= 0.4604*(yup-ylo)+ylo     !   0 sigma,  0.0%
       ysto(8,2,2)= 0.5649*(yup-ylo)+ylo     !   1 sigma, 68.3%
       ysto(9,2,2)= 0.6482*(yup-ylo)+ylo     !                90.0%
c       ysto(10,2,2)= 0.7092*(yup-ylo)+ylo   !   2 sigma, 95.4%
       ysto(10,2,2)= 0.8317*(yup-ylo)+ylo   !   2 sigma, 95.4% sym
       ysto(11,2,2)=0.8489*(yup-ylo)+ylo    !                99.0%
c       ysto(12,2,2)=0.9969*(yup-ylo)+ylo      !   3 sigma, 99.73%
       ysto(12,2,2)=1.2793*(yup-ylo)+ylo      !   3 sigma, 99.73% sym
       ysto(13,2,2)=1.6343*(yup-ylo)+ylo      !   4 sigma, 99.994%

       zlo=out(1,5,i)
       zup=out(2,5,i)

       zsto(1,2,1)=max(0.,0.08764*(zup-zlo)+zlo)   !  -4 sigma, 99.994%
       zsto(2,2,1)= max(0.,0.1779*(zup-zlo)+zlo)   !  -3 sigma, 99.73%
       zsto(3,2,1)= max(0.,0.2128*(zup-zlo)+zlo)   !                99.0%
       zsto(4,2,1)= max(0.,0.2596*(zup-zlo)+zlo)     !  -2 sigma, 95.4%
       zsto(5,2,1)= max(0.,0.2880*(zup-zlo)+zlo)     !                90.0%
       zsto(6,2,1)= max(0.,0.3407*(zup-zlo)+zlo)     !  -1 sigma, 68.3%
       zsto(7,2,1)= max(0.,0.4294*(zup-zlo)+zlo)     !   0 sigma,  0.0%
       zsto(8,2,1)= max(0.,0.5344*(zup-zlo)+zlo)     !   1 sigma, 68.3%
       zsto(9,2,1)= max(0.,0.6156*(zup-zlo)+zlo)     !                90.0%
       zsto(10,2,1)=max(0.,0.6661*(zup-zlo)+zlo)    !   2 sigma, 95.4%
       zsto(11,2,1)=max(0.,0.7597*(zup-zlo)+zlo)      !                99.0%
       zsto(12,2,1)=max(0.,0.8375*(zup-zlo)+zlo)      !   3 sigma, 99.73%
       zsto(13,2,1)=max(0.,1.0656*(zup-zlo)+zlo)      !   4 sigma, 99.994%

       zlo=out(1,6,i)
       zup=out(2,6,i)

       zsto(1,2,2)=max(0.,0.08764*(zup-zlo)+zlo)   !  -4 sigma, 99.994%
       zsto(2,2,2)= max(0.,0.1779*(zup-zlo)+zlo)   !  -3 sigma, 99.73%
       zsto(3,2,2)= max(0.,0.2128*(zup-zlo)+zlo)   !            99.0%
       zsto(4,2,2)= max(0.,0.2596*(zup-zlo)+zlo)   !  -2 sigma, 95.4%
       zsto(5,2,2)= max(0.,0.2880*(zup-zlo)+zlo)   !            90.0%
       zsto(6,2,2)= max(0.,0.3407*(zup-zlo)+zlo)   !  -1 sigma, 68.3%
       zsto(7,2,2)= max(0.,0.4294*(zup-zlo)+zlo)   !   0 sigma,  0.0%
       zsto(8,2,2)= max(0.,0.5344*(zup-zlo)+zlo)   !   1 sigma, 68.3%
       zsto(9,2,2)= max(0.,0.6156*(zup-zlo)+zlo)   !            90.0%
       zsto(10,2,2)=max(0.,0.6661*(zup-zlo)+zlo)   !   2 sigma, 95.4%
       zsto(11,2,2)=max(0.,0.7597*(zup-zlo)+zlo)   !            99.0%
       zsto(12,2,2)=max(0.,0.8375*(zup-zlo)+zlo)   !   3 sigma, 99.73%
       zsto(13,2,2)=max(0.,1.0656*(zup-zlo)+zlo)   !   4 sigma, 99.994%

       gsto(1,2)=max(0.,-0.4812*(gup-glo)+glo)     ! -4 sigma,  99.994%
       gsto(2,2)=max(0.,-0.1542*(gup-glo)+glo)     ! -3 sigma,  99.73%
       gsto(3,2)= max(0.,-0.0406*(gup-glo)+glo)    !            99.0%
       gsto(4,2)= max(0.,0.09895*(gup-glo)+glo)    !  -2 sigma, 95.4%
       gsto(5,2)= max(0.,0.1763*(gup-glo)+glo)     !            90.0%
       gsto(6,2)= max(0.,0.3048*(gup-glo)+glo)     !  -1 sigma, 68.3%
       gsto(7,2)= max(0.,0.4850*(gup-glo)+glo)     !   0 sigma,  0.0%
       gsto(8,2)= max(0.,0.6584*(gup-glo)+glo)     !   1 sigma, 68.3%
       gsto(9,2)= max(0.,0.7752*(gup-glo)+glo)     !            90.0%
       gsto(10,2)=max(0.,0.8433*(gup-glo)+glo)     !   2 sigma, 95.4%
       gsto(11,2)=max(0.,0.9632*(gup-glo)+glo)     !            99.0%
       gsto(12,2)=max(0.,1.0589*(gup-glo)+glo)     !   3 sigma, 99.73%
       gsto(13,2)=max(0.,1.3278*(gup-glo)+glo)     !   4 sigma, 99.994%

c     crossfade from OGE to IGE

       if(out(1,1,i).ge.Tg+1.0) then
        do j=1,13
	 do l=1,2
          ysto(j,1,l)=ysto(j,2,l)
          zsto(j,1,l)=zsto(j,2,l)
          gsto(j,1)=gsto(j,2)
	 enddo
        enddo
       endif

       if(out(1,1,i).gt.Tg-1.0.and.out(1,1,i).lt.Tg+1.0) then
        do j=1,13
	 do l=1,2
          ysto(j,1,l)=(Tg+1.0-out(1,1,i))/2.*ysto(j,1,l)
     &              +(out(1,1,i)-(Tg-1.0))/2.*ysto(j,2,l)
          zsto(j,1,l)=(Tg+1.0-out(1,1,i))/2.*zsto(j,1,l)
     &              +(out(1,1,i)-(Tg-1.0))/2.*zsto(j,2,l)
          gsto(j,1)=(Tg+1.0-out(1,1,i))/2.*gsto(j,1)
     &              +(out(1,1,i)-(Tg-1.0))/2.*gsto(j,2)
         enddo
	enddo
       endif

       prob(1)=99.994
       prob(2)=99.73
       prob(3)=99.0
       prob(4)=95.4
       prob(5)=90.0
       prob(6)=68.3
       prob(7)=0.0
       prob(8)=68.3
       prob(9)=90.0
       prob(10)=95.4
       prob(11)=99.0
       prob(12)=99.73
       prob(13)=99.994

c     for gnuplot 2-sigma, 3-sigma levels

      if (normalize.eq.0) then
         write(20,61)out(1,1,i)*t0,min(ysto(4,1,1)*b0,ysto(4,1,2)*b0),
     &   max(ysto(10,1,1)*b0,ysto(10,1,2)*b0),min(ysto(2,1,1)*b0,
     &   ysto(2,1,2)*b0),max(ysto(12,1,1)*b0,ysto(12,1,2)*b0)
         write(21,61)out(1,1,i)*t0,min(zsto(4,1,1)*b0,zsto(4,1,2)*b0),
     &   max(zsto(10,1,1)*b0,zsto(10,1,2)*b0),min(zsto(2,1,1)*b0,
     &   zsto(2,1,2)*b0),max(zsto(12,1,1)*b0,zsto(12,1,2)*b0)
         write(22,61)out(1,1,i)*t0,gsto(4,1)*gam0,gsto(10,1)*gam0,
     &   gsto(2,1)*gam0,gsto(12,1)*gam0
      else
         write(20,61)out(1,1,i),min(ysto(4,1,1),ysto(4,1,2)),
     &   max(ysto(10,1,1),ysto(10,1,2)),min(ysto(2,1,1),ysto(2,1,2)),
     &   max(ysto(12,1,1),ysto(12,1,2))
         write(21,61)out(1,1,i),min(zsto(4,1,1),zsto(4,1,2)),
     &   max(zsto(10,1,1),zsto(10,1,2)),min(zsto(2,1,1),zsto(2,1,2)),
     &   max(zsto(12,1,1),zsto(12,1,2))
         write(22,61)out(1,1,i),gsto(4,1),gsto(10,1),gsto(2,1),
     &	 gsto(12,1)
      endif	 

c     for gnuplot 0-sigma level

!        write(23,61)out(1,1,i),ysto(7,1,1),ysto(7,1,2),
!      &             zsto(7,1,1),zsto(7,1,2),gsto(7,1)

      enddo

200   continue

      close (10)
!       close (11)
!       close (12)
      close (20)
      close (21)
      close (22)
!       close (23)

c     write met data

      open(11,file='met.dat')

      do k=1,nz
        dtemp=thref-th(k)
        N=0.
        if(dtemp/(z0-zm(k)).gt.0.0) then
          N=(9.81/thref*dtemp/((z0-zm(k))*b0))**0.5*t0
c	N=0.09
        endif
        write(11,61)zm(k),um(k),vm(k),wm(k),qq(k),N,EDRpr(k),th(k)  ! Normalized Qu.
      enddo

      close (11)

c     Scoring
c     Attention: lidar data not normalized !!!!!!!!!!!!!!!
c     gaml(k)/gam0*(A-exp(-C/(nu1*(tfirst/t0-T1))))*gamfac to be used !!!!!!!!!!!!!
c     do not use lidar data with entries = 999 or NaN

       ny = 0
       nz = 0
       nc = 0
       rmsy = 0.
       rmsz = 0.
       rmsc = 0.
       biasyr = 0.
       biaszr = 0.
       biascr = 0.
       biasyl = 0.
       biaszl = 0.
       biascl = 0.



      scordat=case
      write(scordat(22:25),'(A4)')'.dat'
c      print *,scordat
c      open(30,file=scordat)

      do k = ifirst,klidl


c------left----------------------

      tdifold=999.
       do j=1,nn
        tdif=abs(tl(k)/t0-out_mean(1,j))
	if(tdif.lt.tdifold) then
	 ntp=j
         tdifold=tdif
	endif
	if(out_mean(1,j).gt.(tl(k)/t0-0.01)) goto 211
       enddo

211   continue

c       print *,ntpr,'ntpr p2p'
c       print *,ntpl,'ntpl p2p'

c     scoring
     
      if(ntp.gt.0) then

      if (-out_mean(2,ntp).lt.0) goto 212  ! scoring stops for gam<0

       if(abs(yl(k)).gt.0) then             ! makes sure that yl is a number
       	 if (tl(k).ge.0) then
            ny = ny +1
            rmsy = (yl(k)-out_mean(3,ntp)*b0)**2+rmsy
c            print *,out_mean(3,ntp)*b0,'erry1 d2p',yl(k)
c	    print *,out_mean(1,ntp)*t0,tl(k)
            biasyl(k) = (-out_mean(3,ntp)*b0+yl(k))
!             write(24,990)tl(k)/t0,
!      &      yl(k)/b0,out(1,3,ntp),out_mean(3,ntp),out(2,3,ntp)
         endif
      endif


       if(zl(k).gt.0) then
       	 if (tl(k).ge.0) then
        nz = nz +1
        rmsz = (zl(k)-out_mean(4,ntp)*b0)**2+rmsz
        biaszl(k) = (-out_mean(4,ntp)*b0+zl(k))
!         write(25,990) tl(k)/t0,
!      &  zl(k)/b0,out(1,5,ntp),out_mean(4,ntp),out(2,5,ntp)
           endif
       endif

       if(abs(gaml(k)).gt.0) then
        if (tl(k).ge.0) then
         nc = nc +1
         rmsc = (gaml(k)+out_mean(2,ntp)*gam0/gamfac/
     &          (A-exp(-C/(nu1*(tfirst-T1))))/gamfac)**2+rmsc
         biascl(k) = (out_mean(2,ntp)*gam0/gamfac/
     &          (A-exp(-C/(nu1*(tfirst-T1))))/gamfac+gaml(k))
!          write(26,990) tl(k)/t0,
!      &   gaml(k)/gam0*(A-exp(-C/(nu1*(tfirst-T1))))*gamfac**2,
!      &   out(1,2,ntp),-out_mean(2,ntp),out(2,2,ntp)
        endif
       endif


c------right---------------------

      tdifold=999.
       do j=1,nn
        tdif=abs(tr(k)/t0-out_mean(1,j))
	if(tdif.lt.tdifold) then
	 ntp=j
         tdifold=tdif
	endif
	if(out_mean(1,j).gt.(tr(k)/t0-0.01)) goto 213
       enddo

213    continue

       if(abs(yr(k)).gt.0) then
       	 if (tr(k).ge.0) then
           ny = ny +1
           rmsy = (yr(k)-out_mean(6,ntp)*b0)**2+rmsy
c            print *,out_mean(6,ntp)*b0,'erry2 d2p',yr(k)
c	    print *,out_mean(6,ntp)*t0,tr(k)
           biasyr(k) = (-out_mean(6,ntp)*b0+yr(k))
!            write(24,990) tr(k)/t0,
!      &     yr(k)/b0,out(1,4,ntp),out_mean(6,ntp),out(2,4,ntp)
           endif
       endif


       if(zr(k).gt.0) then
       	 if (tr(k).ge.0) then
        nz = nz +1
        rmsz = (zr(k)-out_mean(7,ntp)*b0)**2+rmsz
        biaszr(k) = (-out_mean(7,ntp)*b0+zr(k))
!         write(25,990) tr(k)/t0,
!      &  zr(k)/b0,out(1,6,ntp),out_mean(7,ntp),out(2,6,ntp)
           endif
       endif


       if(abs(gamr(k)).gt.0) then
       	 if (tr(k).ge.0) then
        nc = nc +1
        rmsc = (gamr(k)-out_mean(5,ntp)*gam0/gamfac/
     &         (A-exp(-C/(nu1*(tfirst-T1))))/gamfac)**2+rmsc
        biascr(k) = (-out_mean(5,ntp)*gam0/gamfac/
     &         (A-exp(-C/(nu1*(tfirst-T1))))/gamfac+gamr(k))
!         write(26,990) tr(k)/t0,
!      &  gamr(k)/gam0*(A-exp(-C/(nu1*(tfirst-T1))))*gamfac**2,
!      &  out(1,2,ntp),out_mean(5,ntp),out(2,2,ntp)
           endif
       endif

c-----left/right-----

       if ((abs(yl(k)).gt.0).and.(abs(yr(k)).gt.0)) then
          biasb(k)= -abs(-out_mean(3,ntp)*b0+out_mean(6,ntp)*b0)+
     &              abs(yl(k)-yr(k))	  
       endif

      endif !ntp > 0

      biasyl(k)=biasyl(k)/b0
      biaszl(k)=biaszl(k)/b0
      biascl(k)=biascl(k)/(gam0/gamfac)
      biasyr(k)=biasyr(k)/b0
      biaszr(k)=biaszr(k)/b0
      biascr(k)=biascr(k)/(gam0/gamfac)
      biasb(k)=biasb(k)/b0
!       if (abs(vm(k)).lt.99999) then
!       write(17,194) case,tl(k)/t0,biasyl(k),biaszl(k),
!      &              gaml(k)/gam0*(A-exp(-C/(nu1*(tfirst-T1))))*gamfac,
!      &              biascl(k),tr(k)/t0,biasyr(k),biaszr(k),
!      &              gamr(k)/gam0*(A-exp(-C/(nu1*(tfirst-T1))))*gamfac,
!      &              biascr(k),gam0,biasb(k),wind_dir_lid,z0,
!      &              vm(k)
!       endif
      enddo

212   continue

      rmsy=sqrt(rmsy/ny)/b0
      rmsz=sqrt(rmsz/nz)/b0
      rmsc=sqrt(rmsc/nc)/(gam0/gamfac)


      write(16,995) case,rmsy,rmsz,rmsc
c     write(*,995) case,rmsy,rmsz,rmsc
      write(*,*) '----------------'
      write(*,'(A7,A7,A7)') 'rmsy','rmsz','rmsc'
      write(*,'(F7.2,F7.2,F7.2)') rmsy,rmsz,rmsc
      write(*,*) '----------------'

      close (16)
      close (17)
      close (18)
!       close (24)
!       close (25)
!       close (26)
      close (30)
      close(15)



194   format(A22,15(F10.3))
195   format(15(f10.3))

61    FORMAT(56F10.3)
62    FORMAT(4(A4,F10.3))
63    FORMAT(6A3,I2,35F8.3)
980   format(A3,10(F8.3))
990   format(13(F10.3))
992   format(13(F10.2))
995   format(A22,8(F9.4))
996   format(I4,8(F9.4))
997   format(A6,A3,5(F8.3))
1010   format(2(A8),5(F8.2))
1111   format(A22,2(E11.3),5(F8.3))


      write(*,*) 'ready'

      goto 994

993   print *,'no valid case ', ac, date, time

      write(case(1:20),'(A20)')'xxxxxxxxxxxxxxxxxxxx'

      close (10)
      close (11)
      close (12)
      close (13)
      close (16)
      close (18)

994   continue

!       output='; set output "'
!       write(output(15:23),'(A9)')'P2P_levl_' ! Norm. Qu.
!       write(output(24:44),'(A21)')case
!       write(output(45:48),'(A4)')'.ps"'
! 
!       write(case(7:7),'(A1)')' '
!       write(case(17:17),'(A1)')' '
!       title='; set title   "'//case
      

c      if (normalize.eq.1) then
c        open(12,file='P2P_2sig-levels_SoRa_col_norm.plo'
c     &         ,access='direct',form='formatted',recl=61)
c
c        write(12,'(A61)',rec=1)
c     &  'set terminal postscript portrait enhanced color solid 14'
c        write(12,'(A50)',rec=2) output
c        write(12,'(A39,A4,F6.1,A2,A1)',rec=3) title,
c     &        'z0 =',z0*b0,' m','"'
c
c        close (11)
c        close (12)
c      else
c        open(12,file='P2P_2sig-levels_SoRa_col_dim.plo'
c     &         ,access='direct',form='formatted',recl=61)
c
c        write(12,'(A61)',rec=1)
c     &  'set terminal postscript portrait enhanced color solid 14'
c        write(12,'(A50)',rec=2) output
c        write(12,'(A39,A4,F6.1,A2,A1)',rec=3) title,
c     &        'z0 =',z0*b0,' m','"'
c
c        close (11)
c        close (12)
c      endif

      stop

998   print *,'failed open '
      stop 'open file'

999   write(*,*)'P2P evaluation stopped (erroneous / missing met data)'

      end PROGRAM P2P_REA


      SUBROUTINE PrintIfVerbose(input_string, flag)
            CHARACTER(LEN=*), INTENT(IN) :: input_string
            LOGICAL, INTENT(IN) :: flag
        
            IF (flag) THEN
                PRINT *, input_string
            END IF
        END SUBROUTINE PrintIfVerbose





