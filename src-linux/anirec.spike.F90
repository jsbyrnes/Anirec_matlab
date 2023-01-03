MODULE anirec_module

USE, intrinsic :: iso_c_binding
SAVE

Integer, parameter :: kint    = C_INT
Integer, parameter :: kreal   = C_DOUBLE

Integer(kint), parameter :: ipulse  = 1
Integer(kint), parameter :: per     = 2
Real(kreal)  , parameter :: t1      = 0.00
Real(kreal)  , parameter :: t2      = 100

real(kreal), parameter   :: pi  = 3.14159265358979
real(kreal), parameter   :: eps = 1.d-6
real(kreal), parameter   :: tol = 1.d-3

!real(kreal), dimension(40000) :: a

real(kreal)                    :: w(3,1001),t(3,3),ttl(3,3),s(3,3),stl(3,3),r(3,3),x(3),y(3)
real(kreal), dimension(1000)   :: z,dz,rho,vp,vs,vp2,vp4,vs2,vss
real(kreal), dimension(1000)   :: xmu,xla,xmu2,xla2,xla4
Complex(kreal)                 :: xnu(6,1001),xl(6,1000),pfac(6,3),u(3,6)
real(kreal)                    :: dat44(8200,3,3)
real(kreal)                    :: qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),fv(6)
integer(kint)                  :: iv(6)
real(kreal)                    :: qi(6,6),xr(6),xi(6),yr(6),yi(6)
integer(kint)                  :: ips(3)
Complex(kreal)                 :: pp(3),u0(3),ee(6,6,1001),e1(6,6),e2(6,6),zla(6)
Complex(kreal)                 :: rt(3,3,1001),tt(3,3,1001),rt0(3,3),trc(3,3)
real(kreal)                    :: cc4(8200),zz4(8200),dat4(8200,3,3),ccc(8200)
Complex(kreal)                 :: resp(3,3,2050)
real(kreal)                    :: frqq(2050)
real(kreal)                    :: roota(1001),rootb(1001),jtrval(1001),kroots(2050)
real(kreal)                    :: ievan(10000)
Complex(kreal)                 :: rtm(6,6,1000)
Complex(kreal)                 :: co(6,1001),ur(3)
real(kreal)                    :: idfct(4,1001),adf(2,1001)

Contains

Subroutine anirec ( theta_in, phi_in, z_in, vp_in, A_in, B_in, vs_in, C_in, rho_in, &     ! Model Inputs
                    nlayers, cc_in, baz_in, dt_in, phase,                        &     ! Input parameters
                    z_comp, r_comp, t_comp                                       )  &     ! output arrays
Bind(C, name="anirec")

! This is a mex implementation of the anirec program, program comments preserved.

!  anirec - program to calculate receiver-function response of
!  a stack of anisotropic layers to a plane wave incident from below
!  CAN BE USED IN GRIDSEARCH OVER CC AND BAZ - SEE COMMENTED LINES FOR LOOP
!  cannibalized from aniprop.f  11/18/95
!  xf77 -o anirec_osc -fast -native -O5 anirec.f /data/d4/park/Plotxy/plotlib /data/d4/park/Ritz/eislib /data/d4/park/Ritz/jlib
!  f77 -o anirec -fast -native -O5 anirec.f /data/d4/park/Ritz/eislib /data/d4/park/Ritz/jlib
!  for hexagonally symmetric media
!  reads fast axis orientation, constants A,B,C,D,E from file animodel
!  calculate quadratic eigenvalue problem based on the Christoffel matrix
!  see appendix of P. Shearer's thesis
!
!  read model, phase velocity of incident wave, P, SV, or SH
!
!  calc the eigenvector decomps for the layers
!  loop over frequency, calc reflection/transmission matrices
!  calc 3-comp transfer fct response at surface
!  find distortion of reference wavelet
!

      Implicit None

      real(kreal) :: returnvalue

      !input arrays
      real(kreal), intent( in), dimension(1000)        :: theta_in, phi_in, z_in, vp_in, &
                                                          A_in, B_in, vs_in, C_in, rho_in

      !input parameters
      integer(kint), intent(in)                        :: nlayers, phase
      real(kreal), intent( in)                         :: cc_in, baz_in, dt_in

      ! output 
      real(kreal), intent(out), dimension(2001)       :: z_comp, r_comp, t_comp

      integer(kint) :: i, k, jf, iprint, nfrq, nl, nlm, nlp, npad, npts, npul, nst, lim, iii
      complex(kreal):: zz, z1, z2, z0
      real(kreal)   :: amn,amx,baseline, cccc, cmin, con, cpmin, csmin, df, dur, fac, frqmax, ick
      real(kreal)   :: om, pdelay, phi, radi, rbar, ren, sdelay, theta, time, vbar, xx, cc, baz, dt

      cc = cc_in
      baz = baz_in
      dt = dt_in

      returnvalue = 15

!  we reduce the condition numbers of matrices by
!  normalizing physical quantities to make them dimensionless
!  Ill use the normal-mode normalizations
!  which are a little peculiar for the crust, but what the hey!
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=ren*radi
      con=rbar*vbar**2
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
!  Notes on the angle conventions for w-hat, the axis of symmetry:
!
!  In the anisotropic reflectivity code, subroutine matget *assumes* a
!  coordinate system in which z is down (anti-vertical), x is the radial
!  direction, and y is anti-transverse.  Therefore, the position angles
!  theta,phi for w-hat are tilt relative to down, and azimuth defined as a
!  rotation from x towards y.  This rotation is CCW if viewed from below,
!  and CW if viewed from above.  Since w-hat and -(w-hat) define the same
!  axis of symmetry, the position angles *also* can be defined as
!  theta=(tilt from vertical) and phi=(rotation from anti-x (anti-radial)
!  towards anti-y (transverse)).  Viewed from above, this phi rotation is
!  CW, and defines the strike of w-hat relative to the arrival azimuth of
!  the wave.

!  In order to compute seismograms for a variety of back-azimuths, we
!  assume that the default is a wave approaching from the north, so that
!  radial is south and transverse is east.  Due to this orientation, the
!  synthetic code interprets the layered model as having w-hat position
!  angles defined as theta=(tilt from vertical) phi=(strike CW from N).
!  For an event at back-azimuth psi (CW from N), routine matget rotates w-hat
!  from geographic coordinates to ray-based coordinates before computing
!  reflectivity matrices.  If a wave arrives at back-azimuth psi, the strike
!  of the axis of symmetry w-hat relative to its arrival azimuth is
!  phi'=phi-psi.  The code performs this rotation with this code in
!  subroutine matget, for w-hat azimuth "az":
!
!        caz=dcosd(az)
!        saz=-dsind(az)   ! sin(-az)
!        do n=1,nlp
!          ww(3)=w(3,n)
!          ww(1)=w(1,n)*caz-w(2,n)*saz
!          ww(2)=w(1,n)*saz+w(2,n)*caz
!  ...
!
!  In this manner, the axes of symmetry of the model, saved in array w(.,.),
!  are never modified.
!  in the driver code, "az" is the variable "baz" for back azimuth
!
!  nfrq must be .le.npad/2
      npad=4096
      nfrq=2048
      dur=npad*dt
      df=1./dur
      frqmax=nfrq*df
      !print *,'dt,df,duration of record,max frequency:',dt,df,dur,frqmax
!      print *,'NOTE: cosine^2 taper will be applied up to fmax'
      !print *,'input model? (space-return: animodel)'
      !read(5,102) name
 !102  !format(a)
      !if(name(1:1).eq.' ') then
      !  open(7,file='animodel',form='formatted')
      !else
      !  open(7,file=name,form='formatted')
      !endif
      !read(7,102) title
      !print *,title
      !read(7,*) nl


      !!! ADDED for mex
      nl = nlayers
      !!!

!!!!initalize some arrays

!  read in theta,phi in degrees - polar coords of fast axis
      nlp=nl+1
      nlm=nl-1
      do i=1,nlp

        z(i)=0
        vp(i)=0
        vp2(i)=0
        vp4(i)=0
        vs(i)=0
        vs2(i)=0
        rho(i)=0
        xmu(i)=0
        xmu2(i)=0
        xla(i)=0
        xla2(i)=0
        xla4(i)=0


        !!!Added to account for theta in being an array, same with phi
        theta=theta_in(i)
        phi=phi_in(i)
        !!!

        !!!removed for mex
        !read(7,*) theta,phi
        w(1,i)=dble(sin(theta * pi / 180.0)*cos(phi* pi / 180.0))
        w(2,i)=dble(sin(theta* pi / 180.0)*sin(phi* pi / 180.0))
        w(3,i)=dble(cos(theta* pi / 180.0))
!        print *,(w(j,i),j=1,3)
!  read depth to ith interface, vp (m/sec), pk-to-pk cos(2th) relative P pert
!  pk-to-pk cos(4th) relative P pert, v_s,  pk-to-pk cos(2th) relative S pert
!  density (kg/m**3)

        !!!removed for mex
        !read(7,*) z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i)

        !!!this block added to input the values from matlab
        z(i)=z_in(i)
        vp(i)=vp_in(i)
        vp2(i)=A_in(i)
        vp4(i)=B_in(i)
        vs(i)=vs_in(i)
        vs2(i)=C_in(i)
        rho(i)=rho_in(i)
        !!!

!  recall that we interpret fractional values of b,c,e
!  as peak-to-peak relative velocity perts.
!  therefore, e=0.02 is 2% pert to mu from slowest to fastest
        xmu(i)=rho(i)*vs(i)**2/con
        xmu2(i)=vs2(i)*xmu(i)
        xla(i)=rho(i)*vp(i)**2/con
        xla2(i)=vp2(i)*xla(i)
        xla4(i)=vp4(i)*xla(i)
        vs(i)=vs(i)/vbar
        vp(i)=vp(i)/vbar
        rho(i)=rho(i)/rbar
        z(i)=z(i)/radi
      end do

      !close(7)
      do i=2,nl
        dz(i)=z(i)-z(i-1)
      end do
      dz(1)=z(1)

returnvalue = dz(2)*radi

!  search for cmin, cmax is halfspace velocity
      cmin=vs(1)
      vss(1)=vs(1)
      do i=2,nlp
        if(cmin.gt.vs(i)) cmin=vs(i)
        vss(i)=vs(i)
      end do
  900 csmin=vs(nlp)*vbar/1000.
      cpmin=vp(nlp)*vbar/1000.
!      print *,'minimum phase velocity for S wave (km/sec)',csmin
!      print *,'minimum phase velocity for P wave (km/sec)',cpmin
!      print *,'either BAZ is determined by the loop count'
!      print *,'or, in interactive mode, BAZ is user-input'
!      print *,'if BAZ=0, the wave is arriving from N'
!      print *,'in anirec.bh? the N component is the -RADIAL'
!      print *,'          and the E component is the TRANSVERSE'
!      print *,'if BAZ=90, the wave is arriving from E'
!      print *,'  in anirec.bh? the N component is the -TRANSVERSE'
!      print *,'and the             E component is the -RADIAL'
!      print *,'if BAZ=180, the wave is arriving from S'
!      print *,'  in anirec.bh? the N component is the RADIAL'
!      print *,'and the             E component is the -TRANSVERSE'
!      print *,'if BAZ=270, the wave is arriving from W'
!      print *,'  in anirec.bh? the N component is the TRANSVERSE'
!      print *,'and the             E component is the RADIAL'
!  version for fully interactive input of baz and phase velocity

      !!! Block removed for mex files
      !print *,'enter phase velocity of incident wave (km/sec)'
      !read(5,*) cc
      !if(cc.le.0.d0) go to 950
      !!!

!  non-dimensionalize cc

      cc=cc*1000./vbar

      !!! removed for mex
      !print *,'enter back azimuth of incident wave (deg CW from N)'b
      !read(5,*) baz
      !!!

!  version for grid search over baz and phase velocity
!  we eschew the do-loop indentation convention here
!      do ipvel=10,50,10
!      do ibaz=0,350,10
!      baz=float(ibaz)
!      print *,'phase velocity & back azimuth of incident wave'
!      print *,ipvel,ibaz
!      cc=ipvel*1000./vbar
!  end version

!  print the organ-pipe mode count for 1Hz
!  the lowest layer (nl+1) is taken as evanescent region.
      sdelay=0.
      pdelay=0.
      do i=1,nl
         sdelay=sdelay+(dz(i)/vs(i))/ren
         pdelay=pdelay+(dz(i)/vp(i))/ren
      end do
!      print *, 'organ-pipe mode count at 1 Hz in stack of layers: S & P'
!      print *, 'These are also the S and P travel time from the Moho'
!      print 104,sdelay,pdelay
      
 ! 104 format(2f10.1)
      
!        pause


!  calc the eigenvector decomps for the layersiu
!  need to identify upgoing P,SV,SH in the halfspace
      call matget(nl,cc,baz, returnvalue)

!  loop over frequency, calc reflection/transmission matrices
!  calc 3-comp transfer fct response at surface
      do jf=1,nfrq
        om=2.d0*pi*jf*df/ren
        frqq(jf)=jf*df
        call respget(nl,om,cc,resp(1,1,jf), returnvalue, z_comp)
      end do
      cccc=cc*vbar/1000.

!  plot responses: P,SV,SH
!      do ity=1,3
!        do k=1,3
!          do i=1,nfrq
!            cc4(i)=zabs(resp(k,ity,i))
!          end do
!          call plotit(frqq,cc4,dum,nfrq,title,xlabel,
!     x      ylabel(ity),2,0,0.05*(3-k),0,(30+ity)*(k/3))
!        end do
!      end do
!  lets run a pulse thru these functions
!  s(t)=cos(2pi*t/T)*sin^2(2pi*t/2T) for two oscillations of the cos
!  OR
!  s(t)=sin(2pi*t/2T)**2 for 1/2 oscillation of the sin
!  OR bbpulse
!  s(t)=t*cos(2.*pi*t/sqrt(1.+t**2))*exp(-0.4*t**2)
!  let T=1 and 2 sec
! version for user-supplied cc and baz
! version for gridsearch cc and baz
!      ipulse=3
!
!  another version for gridsearch
!      ipulse=2
!      per=8.
! end gridsearch lines
      dur=npad*dt
      !t1=0.
      !t2=dur
! version for usersupplied cc and baz
! version for gridsearch cc and baz
!      t1=0.
!      t2=60.
!  end versions
      npts=t2/dt+1
      nst=t1/dt+1
!  we start at dt, with duration 2T -- ONE CYCLE
      npul=per/dt
      do i=1,npad
        zz4(i)=0.d0
      end do
      fac=2.d0*pi/per
      if(ipulse.eq.1) then
        do i=1,npul
          time=i*dt
          zz4(i)=(dsin(fac*time/2.d0))**2
        end do
      elseif(ipulse.eq.2) then
        do i=1,npul
          time=i*dt
          zz4(i)=dsin(fac*time)*(dsin(fac*time/2.d0))**2
        end do
      elseif(ipulse.eq.3) then
        do i=1,160
          xx=0.05*(i-80)
          zz4(i)=xx*cos(2.*pi*xx/sqrt(1.+xx**2))*exp(-0.4*xx**2)
        end do
      else
        do i=1,160
          xx=0.05*(i-80)
          zz4(i)=cos(5.*pi*xx/sqrt(1.+4*xx**2))*exp(-4.0*xx**2)
        end do
      endif
!      print *,(zz4(i),i=1,npul)

call refft(zz4,npad,1,1)

!  zero the DC and Nyquist
!  ick switches the sign of y and z components to transverse & vertical
      zz4(1)=0.
      zz4(2)=0.
      do iprint=phase, phase
        ick=1
        baseline=0.
        do k=1,3
          cc4(1)=0.
          cc4(2)=0.
          if(k.gt.1) ick=-1
          do jf=1,nfrq
!           zz=ick*dcmplx(dble(zz4(2*jf+1)),dble(zz4(2*jf+2)))
	       zz=ick
            zz=zz*resp(k,iprint,jf)
            cc4(2*jf+1)=dreal(zz)
            cc4(2*jf+2)=dimag(zz)
          end do
          do jf=2*nfrq+3,npad
            cc4(jf)=0.
          end do
          call refft(cc4,npad,-1,-1)
          amx=cc4(1)
          amn=cc4(1)
          do i=1,npts
            amx=amax1(amx,cc4(i+nst))  
            amn=amin1(amn,cc4(i+nst))  
          end do
          baseline=baseline-amn
          do i=1,npad
            dat4(i,k,iprint)=cc4(i)+baseline
            dat44(i,k,iprint)=cc4(i)
          end do
          baseline=baseline+amx
        end do
      end do

lim = t2/dt+1

do i=1,lim

r_comp(i)=dat44(i, 1, phase)
!r_comp(i)=cc4(i)
t_comp(i)=dat44(i, 2, phase)
z_comp(i)=dat44(i, 3, phase)
end do

returnvalue = 55


      !!!
      !Don't want to use the old sac routines, I want to save the traces in 
      !in z_comp, r_comp, t_comp
      !!!

!      cc4(1)=nst*dt
!      cc4(2)=dt
      
!      wf(1) = '.P'
!      wf(2) = '.SV'
!      wf(3) = '.SH'
      
      !print *, sdelay, pdelay, cccc
!      do ity=1,1
          !name = 'split.0.0.bhz'
          !call sacread(name,a,ierr)
          !if(ierr.ne.0) go to 950    ! in case of misprintd filename
          !ahead(1)=dt                !change dt in sac header

! use organ-pipe count to correct for traveltime of Swave thru stack
!          tdelay=sdelay-13.
! use organ-pipe count to correct for traveltime of Pwave thru stack
!          tdelay=pdelay-4.
	  
! note 13 and 4 are the time shift from zero for S and P.
!	    if (ity .ne.1) then
!	     tdelay =sdelay - 10.
!	    else
!	     tdelay=pdelay - 10.
!	    endif
!          nst0=tdelay/dt
! let's find the real zero time and shift to +10 sec
          
!	    amx = 0.
!            do i = 2, npad
!	       if (ity .eq. 1) amn=dat44(i,3,1)
!	       if (ity .eq. 2) amn=dat44(i,1,2)
!	       if (ity .eq. 3) amn=dat44(i,2,3)
		   
!	       if (amx .lt. abs(amn) ) then
!	           amx = amn
!	           nst0 = i
!	           nst0 = i
!	       endif
!         enddo
 !         print *, nst0
!          shift zero to +10 sec
!          if(nst0.lt.0) nst0=0
!          nst0 = nst0 - 5./dt
! shift data to : 5 sec for P incident wave; center for S waves
!	        if ((nst0 - npts/2) .lt.0) then
!	           ahead(47) = nst0 *dt
!			   nst0 =0
!              else
!	           ahead(47) = npts/2 *dt
!			   nst0 = nst0 - npts/2
!	        endif


! comp1=radial comp2=transverse comp3=vert
!  we transform radial & transverse into east & north
!          ahead(53)=baz
!	    ahead(46)=1./(cc/1000.*vbar)
	  
!          ahead(53)=baz
!  here we rotate the horizontals
!  if we are in the loop, BAZ is determined by the loop count
!  if we are in interactive mode, BAZ=0, and the wave is arriving from N
!  if BAZ=0, the wave is arriving from N
!    in anirec.bh? the N component is the -RADIAL
!              and the E component is the TRANSVERSE
!  if BAZ=90, the wave is arriving from E
!    in anirec.bh? the N component is the -TRANSVERSE
!  and the             E component is the -RADIAL
!  if BAZ=180, the wave is arriving from S
!    in anirec.bh? the N component is the RADIAL
!  and the             E component is the -TRANSVERSE
!  if BAZ=270, the wave is arriving from W
!    in anirec.bh? the N component is the TRANSVERSE
!  and the             E component is the RADIAL
!          cs=cos(baz* pi / 180.0)
!          sn=sin(baz* pi / 180.0)
!          iah(80) = npts        !number of points

!          iseed=1234
!          who=rand(iseed)

!          write(name2,'(a,f5.2,a,i3,a,a)') 'anif.',
!              cc/1000*vbar,'.',int(baz),'.000',wf(ity)
!          do iii=1,18
!             if (name2(iii:iii).eq.' ') name2(iii:iii)='0'
!          end do
	  
!          name2(15:18)='.bhr'

!  write BOTH radial/transverse and north/east horizontal components
!  the radial/transverse SAC files are less vulnerable to BAZ confusion in SAC


!          ahead(58)=baz+180.
!          ahead(59)=90.
!          chead(151)='BHR '
!	  do i=1,nst
!            a(i)=0.0*rand(0)

!             a(i) = 0
!	  end do
!          do i=nst,iah(80)
 !           a(i)=dat44(nst0+i-nst,1,ity)
!            a(i)=dat44(nst0+i-nst,1,ity)+0.0*rand(0)
!          end do
!          call sacout(name2,a)

!          name2(15:18)='.bht'
!          ahead(58)=baz+90.
!          ahead(59)=90.
!          chead(151)='BHT '
!	      do i=1,nst
!               a(i)=0.00*rand(0)
!	      enddo
!          do i=nst,iah(80)
!            a(i)=dat44(nst0+i-nst,2,ity)+0.0*rand(0)
!          end do
 	  
!c          do i=1,iah(80)
!c            a(i)=dat44(nst0+i,2,ity)+0.01*rand(0)
!c          end do
!          call sacout(name2,a)

!          name2(15:18)='.bhz'
!          ahead(58)=0.
!          ahead(59)=0.
!          chead(151)='BHZ '

!          do i=1,nst
!            a(i)=0.00*rand(0)
!	      end do
!          do i=nst,iah(80)
!            a(i)=dat44(nst0+i-nst,3,ity)+0.0*rand(0)
!          end do
!
!c          do i=1,iah(80)
!c            a(i)=dat44(nst0+i,3,ity)+0.01*rand(0)
!c          end do
!          call sacout(name2,a)

!c          name(48:48)='n'
!c          write(name,103) 'anirec.bhn'
!c          ahead(58)=0.
!cc          ahead(59)=90.
!c          chead(151)='BHN '
!c          do i=1,iah(80)
!c            a(i)=-sn*dat44(nst0+i,2,ity)-cs*dat44(nst0+i,1,ity)
!c c    x              					+0.01*rand(0)
!c          end do
!c          print 101, (chead(i),i=111,158)
!c          call sacout(name2,a)
!c          name(48:48)='e'
!c          write(name,103) 'anirec.bhe'
!c          ahead(58)=90.
!c          ahead(59)=90.
!c          chead(151)='BHE '
!c          do i=1,iah(80)
!c            a(i)=cs*dat44(nst0+i,2,ity)-sn*dat44(nst0+i,1,ity)
!c     x              					+0.01*rand(0)
!c          end do
!c          print 101, (chead(i),i=111,158)
!c          call sacout(name,a)

	
!     end do
!  103 format(a,i2,a,i3,a)
!c version for user-supplied cc and baz
!      go to 960
!c end version for user-supplied cc and baz
!c version for grid search over cc and baz
!c      end do
!c      end do
!c  end version for grid search over cc and baz
!  950 continue
!  101 format(80a)
!      stop
end

subroutine matget(nl,cc,az, returnvalue)

!c  SPECIAL VERSION: az rotates the w-hat vector by -az degrees
!c  returns stress-displacement vectors for a stack of anisotropic layers
!c  P waves may be evanescent, but S waves are oscillatory in the stack
!c  the weirdness seen in the surface wave code should not appear in
!c  a receiver function code
!c  however, the iev parameter is retained to avoid leaving timebombs

      implicit none

      complex(kreal)            :: pw,uw,pu
      complex(kreal)            :: z1, z2, eye, z0

      real(kreal)               :: ww(3)

      integer(kint)             :: nl, i, iev, iflag, j, k, jdown, jup, ki, kk, n, nlm, nn, nlp
      integer(kint)             :: is1, is2, iii, jjj, kkk, ierr
      real(kreal)               :: cc, a, b, c, d, e, abcde, bce, az, caz, con, de, fac, facs
      real(kreal)               :: facr, fact, isum, px, radi, rbar, ren, saz
      real(kreal)               :: sum, test, vbar, wii, wr0, wrr

real(kreal), intent(inout) :: returnvalue

!c  set iev=1 ** should be superfluous, save for now
!c  toggle to iev=0 if there is a purely propagating wave in the top layer n=1
      iev=1
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=radi*ren
      con=rbar*radi*radi*ren*ren
      nlp=nl+1
      nlm=nl-1
!c  first calculate vertical wavenumbers and propagating waves for each layer
!c   requires an eigenvector problem be solved
!c  in general, the evanescent vertical wavenumbers have nonzero real parts
!c   complex exponential fct is used to avoid endless branching
!c  horizontal slowness p_x
      px=1.d0/cc
      caz=dble(cos(az* pi / 180.0))
      saz=-dble(sin(az* pi / 180.0))   ! sin(-az)

!!!!!!!!
!changed for mex, loop is infinite if you do over the limit
!!!!!!!!
      do n=1,nlp
        ww(3)=w(3,n)
        ww(1)=w(1,n)*caz-w(2,n)*saz
        ww(2)=w(1,n)*saz+w(2,n)*caz
        a=xla(n)
        b=xla2(n)
        c=xla4(n)
        d=xmu(n)
        e=xmu2(n)

!c        print *,'a,b,c,d,e',a,b,c,d,e
        fact=8.d0*ww(1)*ww(1)*c+2.d0*e
        facs=16.d0*ww(1)*ww(3)*c
        facr=8.d0*ww(3)*ww(3)*c+2.d0*e
!c        print *,'a,b,c,d,e',a,b,c,d,e
!c        print *,'w(.,n),fact,facs,facr',(ww(l),l=1,3),fact,facs,facr
        do i=1,3
!c  first the what-0-what tensor
          do j=1,3
            t(j,i)=fact*ww(j)*ww(i)
            s(j,i)=facs*ww(j)*ww(i)
            r(j,i)=facr*ww(j)*ww(i)
          end do
!c  next the identity tensor - correct an error on 7/6/95
          t(i,i)=t(i,i)+d+e*(2.d0*ww(1)*ww(1)-1.d0)
          s(i,i)=s(i,i)+4.d0*e*ww(1)*ww(3)
          r(i,i)=r(i,i)+d+e*(2.d0*ww(3)*ww(3)-1.d0)
        end do

!c        print 101,(ww(i),i=1,3)
!c        print 101,fact,facs,facr
!c        print *,'t,s,r'
!c        print 101,((t(i,j),j=1,3),i=1,3)
!c        print 101,((s(i,j),j=1,3),i=1,3)
!c        print 101,((r(i,j),j=1,3),i=1,3)
        fac=b-4.d0*c-2.d0*e
!c  next the what-0-xhat and what-0-zhat tensors
        do i=1,3
          t(1,i)=t(1,i)+fac*ww(1)*ww(i)
          t(i,1)=t(i,1)+fac*ww(1)*ww(i)
          s(1,i)=s(1,i)+fac*ww(3)*ww(i)
          s(i,1)=s(i,1)+fac*ww(3)*ww(i)
          s(3,i)=s(3,i)+fac*ww(1)*ww(i)
          s(i,3)=s(i,3)+fac*ww(1)*ww(i)
          r(3,i)=r(3,i)+fac*ww(3)*ww(i)
          r(i,3)=r(i,3)+fac*ww(3)*ww(i)
        end do

        fac=a-b+c-d+e

!c  finally the xhat-0-xhat, zhat-0-zhat, xhat-0-zhat, zhat-0-xhat tensors
        t(1,1)=t(1,1)+fac
        s(3,1)=s(3,1)+fac
        s(1,3)=s(1,3)+fac
        r(3,3)=r(3,3)+fac
!c  mult by horizontal slowness and calc the modified T-matrix
        do i=1,3
          do j=1,3
            t(j,i)=t(j,i)*px*px
            s(j,i)=s(j,i)*px
          end do
          t(i,i)=t(i,i)-rho(n)
        end do

!c  calculate R**(-1).S, R**(-1).T, using routine solve

        nn=3
        do i=1,3
          do j=1,3
            y(j)=s(j,i)
          end do
          call solve(nn,r,x,y, returnvalue)
          do j=1,3
            stl(j,i)=x(j)
          end do
          nn=-3
        end do

        do i=1,3
          do j=1,3
            y(j)=t(j,i)
          end do
          call solve(nn,r,x,y, returnvalue)
          do j=1,3
            ttl(j,i)=x(j)
          end do
        end do
!c  fill the 6x6 Q-matrix
        do i=1,3
          do j=1,3
            qq(j,i)=-stl(j,i)
            qq(j,i+3)=-ttl(j,i)
            qq(j+3,i)=0.d0
            qq(j+3,i+3)=0.d0
          end do
          qq(i+3,i)=1.d0
        end do
!c  solve eigenvalue problem for polarization vectors and vertical slownesses
!c  matrix system is nonsymmetric real valued
!c  solution from the eispack guide

        call balanc(6,6,qq,is1,is2,fv)
        call elmhes(6,6,is1,is2,qq,iv)
        call eltran(6,6,is1,is2,qq,iv,zr)
        call hqr2(6,6,is1,is2,qq,wr,wi,zr,ierr)
!        if(ierr.ne.0) then
!          print *, ierr,'   error!'
!          stop
!        endif
        call balbak(6,6,is1,is2,fv,6,zr)

!c        print *,'for layer',n
!c        print *, 'for phase velocity',cc,'  the vertical slownesses are'
!c        print 101,(wr(i),wi(i),i=1,6)
!c        pause
!  101 format(6g12.4)
!c  eigenvector unpacking, see EISPACK guide, page 88
!c  bad eigenvector order is flagged by wi(i)>0. for odd i
        iflag=0
        do i=1,6
          if(wi(i).eq.0.d0) then
            if(n.eq.1) iev=0
            do j=1,6
              zi(j,i)=0.d0
            end do
          elseif(wi(i).gt.0.d0) then
!c  bad eigenvector order is flagged by wi(i)>0 for even i
            if((i/2)*2.eq.i) then
              iflag=iflag+1
              iv(iflag)=i
            endif
            do j=1,6
              zi(j,i)=zr(j,i+1)
            end do
          else
            do j=1,6
              zi(j,i)=-zi(j,i-1)
              zr(j,i)=zr(j,i-1)
            end do
          endif
!c  normalize by the last three indices
          sum=0.d0
          do j=4,6
            sum=sum+zr(j,i)**2+zi(j,i)**2
          end do
          sum=dsqrt(sum)
          do j=1,6
            zr(j,i)=zr(j,i)/sum
            zi(j,i)=zi(j,i)/sum
          end do         
        end do
!c  assemble the stress-displacement vectors
!c  calculate the traction components, with i removed
        pp(1)=dcmplx(px,0.d0)
        pp(2)=z0
        do k=1,6
          pp(3)=dcmplx(wr(k),wi(k))
          do i=1,3
            u0(i)=dcmplx(zr(i+3,k),zi(i+3,k))
          end do
          pu=z0
          pw=z0
          uw=z0
          abcde=a-b+c-2.d0*d+2.d0*e
          bce=b-4.d0*c-4.d0*e
          de=d-e
          do i=1,3
            pu=pu+pp(i)*u0(i)
            pw=pw+pp(i)*ww(i)
            uw=uw+u0(i)*ww(i)
          end do
          do i=1,3
            e1(i,k)=u0(i)
            e1(i+3,k)=ww(i)*(pu*ww(3)*bce+8.d0*pw*uw*ww(3)*c+2.d0*(pw*u0(3)+uw*pp(3))*e)
            e1(i+3,k)=e1(i+3,k)+pp(i)*(u0(3)*de+2.d0*uw*ww(3)*e)
            e1(i+3,k)=e1(i+3,k)+u0(i)*(pp(3)*de+2.d0*pw*ww(3)*e)
          end do
          e1(6,k)=e1(6,k)+pu*abcde+pw*uw*bce
!c  almost lastly, mult traction by i
          do i=1,3
            e1(i+3,k)=eye*e1(i+3,k)
          end do
        end do
!c  reorder into upgoing and downgoing waves
!c  we use the exp(-i*omega*t) convention with z increasing downward
!c  so downgoing oscillatory waves have p_z>0, k_z real
!c  downgoing evanescent waves have Im(p_z)>0
!c  if the axis of symmetry is tilted, there are cases where a pair of
!c  near-horizontal plane waves will be both upgoing or both downgoing
!c  since Chen's algorithm depends on a 3,3 split, we must adopt a kluge
!c  similarly, there are cases where the EISPACK routines dont return
!c  the vertical wavenumbers in ordered pairs, but mix them up a bit
!c  this seems to cause problems, so a fix is necessary
!c
!c  first, test for bad eigenvector order, switch k-1->k+1, k->k-1, k+1->k
!c   worst case is iflag=2, real,imag1+,imag1-,imag2+,imag2-,real
        if(iflag.gt.0) then
          do i=1,iflag
            k=iv(i)
            wrr=wr(k-1)
            wii=wi(k-1)
            wr(k-1)=wr(k)
            wi(k-1)=wi(k)
            wr(k)=wr(k+1)
            wi(k)=wi(k+1)
            wr(k+1)=wrr
            wi(k+1)=wii
            do j=1,6
              pu=e1(j,k-1)
              e1(j,k-1)=e1(j,k)
              e1(j,k)=e1(j,k+1)
              e1(j,k+1)=pu
            end do
          end do
        endif
!c  second, divide into upgoing and downgoing waves
        isum=0
        do k=1,6
          iv(k)=0
          if(wi(k).eq.0.d0.and.wr(k).gt.0) iv(k)=1
          if(wi(k).gt.0.d0) iv(k)=1
          isum=isum+iv(k)
        end do
!c  if up and downgoing cohorts are not equal, switch the sense of the
!c  pure-oscillatory wave with smallest wavenumber
  140   continue
        if(isum.ne.3) then
          wr0=0.d0
          do k=1,6
            wr0=dmax1(wr0,dabs(wr(k)))
          end do
          do k=1,6
            if(wi(k).eq.0.d0) then
              if(dabs(wr(k)).lt.wr0) then
                wr0=dabs(wr(k))
                kk=k
              endif
            endif
          end do
          if(iv(kk).eq.0) then
            iv(kk)=1
          else
            iv(kk)=0
          endif
!c  check that we have equal up/down cohorts
          isum=0
          do k=1,6
            isum=isum+iv(k)
          end do
          go to 140
        endif
        jdown=1
        jup=4
!c        print *,'for layer',n,'  the vert wavenums are (0=up,1=dn)'
        do k=1,6
          if(iv(k).eq.1) then
            ki=jdown
            jdown=jdown+1
          else
            ki=jup
            jup=jup+1
          endif
          do i=1,6
            ee(i,ki,n)=e1(i,k)
          end do
!c  incorporate the factor of i into the stored vertical slowness
          xnu(ki,n)=dcmplx(-wi(k),wr(k))
        end do

! 1008 format(a,2g15.6,a,2g15.6)
      end do
!c  now, must identify which upgoing waves in the halfspace are P,SV,SH
!c  crud, this goes back to array ee
!c  3: SH is y-motion
!c  2: SV is (-sqrt((1/vs)**2-p_x**2),0,-p_x)  ! recall that z points down
!c  1: P is (p_x,0,-sqrt((1/vp)**2-p_x**2)
!c  so we branch on size of u_y, and relative sign of u_x and u_z
!c      print *,'in the halfspace:'
!c      do i=4,6
!c        print *,'for i*k_z=',xnu(i,nlp),', the disp-stress vector is'
!c        do j=1,6
!c          xi(j)=dimag(ee(j,i,nlp))
!c          xr(j)=dreal(ee(j,i,nlp))
!c        end do
!        print 101,(xr(j),j=1,6),(xi(j),j=1,6)
!      end do
      do i=4,6
        ips(i-3)=3
        if(zabs(ee(2,i,nlp)).lt.dsqrt(tol)) then  ! not SH
          test=dreal(ee(1,i,nlp))/dreal(ee(3,i,nlp))
          if(test.gt.0.d0) then
            ips(i-3)=2
          else
            ips(i-3)=1
          endif
        endif
      end do
!      print *,'wave prints:',(ips(i),i=1,3)
      return

end

subroutine respget(nl,om,cc,resp2, returnvalue, z_comp)

!c  returns surface response for a stack of anisotropic layers
!c  incident p,sv,sh waves with freq om and phase velocity cc
!c  iev=1 if the waves are evanescent in the top layer, iev=0 otherwise

      implicit none

      integer(kint)                 :: nl, i, iev, ii, iup, j, jj, k, k1, n, nlm, nlp, nn, iii, jjj, kkk
      integer(kint)                 :: k2
      real(kreal)                   :: om, cc, con, h1, h2, radi, ren, vbar, rbar
      complex(kreal)                :: z1,z0,eye
    !!! resp renamed resp2 while local, so it doesn't conflict
      complex(kreal), intent(inout) :: resp2(3,3)

real(kreal), intent(inout) :: returnvalue
real(kreal), intent(inout) :: z_comp(*)

!  set iev=1
!c  toggle to iev=0 if there is a purely propagating wave in the top layer n=1
      iev=1
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=radi*ren
      con=rbar*radi*radi*ren*ren
      nlp=nl+1
      nlm=nl-1
!c  first calculate vertical wavenumbers and propagating waves for each layer
!c    an eigenvector problem was solved in prior subroutine
!c   results stored in array ee(6,6,101)
!c  in general, the evanescent vertical wavenumbers have nonzero real parts
!c   complex exponential fct is used to avoid endless branching
!c
!c  calculate modified R/T coefficients
!c  first calc the propagation factors
!c  note that for dipping fast axes the upgoing and downgoing wavenumbers are
!c  independent, so we must calc all to be safe
      do n=1,nl
        do k=1,3
          xl(k,n)=zexp(om*xnu(k,n)*dz(n))		! downgoing
          xl(k+3,n)=zexp(-om*xnu(k+3,n)*dz(n))	! upgoing
        end do
      end do

!c      do i=1,6
!c        print 1002,xnu(i,3),xl(i,3)
!c      end do
! 1002 format('i*k_z:',2g15.6,',  propfac is',2g15.6)
!c  calculate modified R/T coefficients at each interface

iii = 0

      do n=1,nl
!c  rearrange to e1: waves approaching and e2: waves leaving an interface
        do k=1,3
          do i=1,6
            e1(i,k)=ee(i,k,n+1)
            e2(i,k)=ee(i,k,n)
            e1(i,k+3)=-ee(i,k+3,n)
            e2(i,k+3)=-ee(i,k+3,n+1)
          end do
          zla(k)=xl(k,n)
          if(n.lt.nl) then
            zla(k+3)=xl(k+3,n+1)
          else

!c  reference the upcoming wave amplitude to the top of halfspace
!c  therefore no propagation factor, not valid for evanescent waves in halfspace
!c  in surface wave code this is zero, so that upgoing evanescent waves vanish
            zla(k+3)=1.d0
          endif
        end do

!c mult the columns of e2
        do k=1,6
          do i=1,6
            e2(i,k)=e2(i,k)*zla(k)
          end do
        end do
!c  the possibility of defective matrices must be contemplated here
!c  k=1,2,3 columns are downgoing in nth layer
!c  k=4,5,6 columns are upgoing in (n+1)th layer
!c  the vector e2(.,k1) has already been multiplied by exponential factor zla
        if(idfct(1,n).ne.0) then
          k1=idfct(1,n)
          k2=idfct(2,n)
          do i=1,6
            e2(i,k2)=e2(i,k2)+adf(1,n)*dz(n)*e2(i,k1)
          end do
        endif
!c  the sign change on dz is for upgoing waves
        if(idfct(3,n+1).ne.0) then
          k1=idfct(3,n+1)
          k2=idfct(4,n+1)
          do i=1,6
            e2(i,k2)=e2(i,k2)-adf(2,n+1)*dz(n+1)*e2(i,k1)
          end do
        endif
!c  in order to use csolve to invert e1, must separate into real/imag parts
!c  its clumsy, but im lazy
!c  we calc e1^{-1}\cdot e2\cdot \Gamma one column at a time
        do k=1,6
          do i=1,6
            qq(i,k)=dreal(e1(i,k))
            qi(i,k)=dimag(e1(i,k))
          end do
        end do
        nn=6
        do k=1,6
          do i=1,6
            yr(i)=dreal(e2(i,k))
            yi(i)=dimag(e2(i,k))
          end do
          call csolve(nn,qq,qi,xr,xi,yr,yi, returnvalue)
          nn=-6
          do i=1,6
            rtm(i,k,n)=dcmplx(xr(i),xi(i))
          end do
        end do  
      end do

!c  calc R_ud at the free surface
!c  note that first two factors in Chen (20) dont collapse
!c  mult by inv-matrix one column at a time
      do k=1,3
        do i=1,3
          rt0(i,k)=ee(i+3,k+3,1)*xl(k+3,1)
          s(i,k)=dreal(ee(i+3,k,1))
          t(i,k)=dimag(ee(i+3,k,1))
        end do
      end do
!c  the possibility of defective matrices must be contemplated here
!c  these waves are upgoing in 1st layer
!c  the sign change on dz is for upgoing waves, and xl(k1,1)=xl(k2,1)
      if(idfct(3,1).ne.0) then
        k1=idfct(3,1)
        k2=idfct(4,1)-3
        do i=1,3
          rt0(i,k2)=rt0(i,k2)-adf(2,1)*dz(1)*ee(i+3,k1,1)*xl(k1,1)
        end do
      endif
      nn=3
      do k=1,3
        do i=1,3
          yr(i)=dreal(rt0(i,k))
          yi(i)=dimag(rt0(i,k))
        end do
        call csolve(nn,s,t,xr,xi,yr,yi, returnvalue)
        nn=-3
        do i=1,3
          rt0(i,k)=-dcmplx(xr(i),xi(i))
        end do
      end do

!c  recursive calc of generalized R/T coefs:
!c  in contrast to the surface-wave code, we start from the top layer and
!c  iterate down to the halfspace
!c  we also uses submatrices of generalized R/T matrix in different order
      do n=1,nl
!c  first the generalized upward-transmission coef:
        do k=1,3
          do i=1,3
            trc(i,k)=z0
            if(n.gt.1) then
              do j=1,3
                trc(i,k)=trc(i,k)-rtm(i+3,j,n)*rt(j,k,n-1)
              end do
            else
!c  use free-surface reflection matrix in top layer (interface "zero")
              do j=1,3
                trc(i,k)=trc(i,k)-rtm(i+3,j,n)*rt0(j,k)
              end do
            endif
          end do
          trc(k,k)=trc(k,k)+z1
        end do  
        do k=1,3
          do i=1,3
            s(i,k)=dreal(trc(i,k))
            t(i,k)=dimag(trc(i,k))
          end do
        end do
        nn=3
        do k=1,3
          do i=1,3
            yr(i)=dreal(rtm(i+3,k+3,n))
            yi(i)=dimag(rtm(i+3,k+3,n))
          end do
          call csolve(nn,s,t,xr,xi,yr,yi, returnvalue)
          nn=-3
          do i=1,3
            tt(i,k,n)=dcmplx(xr(i),xi(i))
          end do
        end do  
!c  next the generalized reflection coef:
        do k=1,3
          do i=1,3
            trc(i,k)=z0
            if(n.gt.1) then
              do j=1,3
                trc(i,k)=trc(i,k)+rt(i,j,n-1)*tt(j,k,n)
              end do
            else
!c  use free-surface reflection matrix in top layer (interface "zero")
              do j=1,3
                trc(i,k)=trc(i,k)+rt0(i,j)*tt(j,k,n)
              end do
            endif
          end do
        end do  
        do k=1,3
          do i=1,3
            rt(i,k,n)=rtm(i,k+3,n)
            do j=1,3
              rt(i,k,n)=rt(i,k,n)+rtm(i,j,n)*trc(j,k)
            end do
          end do
        end do  
      end do

! 1001 format(6f14.6)
!c      print *,'free-surface reflection'
!c      print 1001,((rt0(i,j),j=1,3),i=1,3)
!c      do n=1,nl
!c        print *,'interface',n
!c        print 1001,((rt(i,j,n),j=1,3),i=1,3)
!c        print 1001,((tt(i,j,n),j=1,3),i=1,3)
!c      end do
!c  using the p,sv,sh identification, we propagate upward to the surface,
!c  calculate the wave coefs in the top layer, then the particle displacement
      do iup=1,3
        do i=1,3
          co(i+3,nlp)=z0
        end do
        co(iup+3,nlp)=z1
!c  from upgoing coefs in the n+1 layer, calculate
!c  upgoing coefs in the nth layer, downgoing coefs in the n+1 layer
        do n=nl,1,-1
          do i=1,3
            co(i+3,n)=z0
            co(i,n+1)=z0
            do j=1,3
              co(i+3,n)=co(i+3,n)+tt(i,j,n)*co(j+3,n+1)
              co(i,n+1)=co(i,n+1)+rt(i,j,n)*co(j+3,n+1)
            end do
          end do
        end do
!c  then downgoing coefs in the top layer:
        do i=1,3
          co(i,1)=z0
          do j=1,3
            co(i,1)=co(i,1)+rt0(i,j)*co(j+3,1)
          end do
        end do
!c        print *,'upgoing coefs'
!c        print 1001,((co(j+3,n),j=1,3),n=1,nlp)
!c        print *,'downgoing coefs'
!c        print 1001,((co(j,n),j=1,3),n=1,nlp)
!c  calc the surface displacement
        h1=0.d0
        h2=z(1)
        do i=1,3
          ur(i)=z0
          do k=1,3
            ur(i)=ur(i)+co(k,1)*ee(i,k,1)+co(k+3,1)*ee(i,k+3,1)*(zexp(om*xnu(k+3,1)*(-h2)))
          end do
!c  check for the xtra terms associated with defective matrices
          if(idfct(1,1).ne.0) then
            ii=idfct(1,1)
            jj=idfct(2,1)
            ur(i)=ur(i)+co(jj,1)*adf(1,1)*ee(i,ii,1)*(-h1)
          endif
          if(idfct(3,1).ne.0) then
            ii=idfct(3,1)
            jj=idfct(4,1)
            ur(i)=ur(i)+co(jj,1)*adf(2,1)*ee(i,ii,1)*(-h2)*(zexp(om*xnu(ii,1)*(-h2)))
          endif
        end do
!c  copy the surface displacement into the response matrix

        do i=1,3
          resp2(i,ips(iup))=ur(i)
        end do

      end do

      return
end

#ifdef NOT_USED

!!! removed adf from decleration
subroutine defective(i,j,n,a,b,c,d,e,px)


!c  kluge for dealing with nearly defective propagator matrices
!c  in which the eigenvectors,
!c  which represent the particle motion of upgoing and downgoing waves
!c  become nearly parallel.
!c  in this case the solution for system of ODEs is
!c  a_1 \bf_1 e^xnu*(z-z0) + a_2*(\bf_2 + adf*(z-z0)*\bf_1)e^xnu*(z-z0)

      implicit none

      integer(kint)            :: i, j, n

      complex(kreal)           :: z1,z0,znu,eye
      complex(kreal)           :: zq1(3,3),zq2(3,3),u1(3),u2(3),zq3(2,2),xee(3)
      real(kreal)              :: q1r(3,3),q1i(3,3),q2r(3,3),q2i(3,3),fv2(3),fv3(3)

      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
!c  for the extravector, need to solve system of equations
!c  based on original 6x6 Q matrix
!c  the plane-wave solutions generalize to the form
!c  u0*e^{i*nu*(z-z0)}  and  u1*e^{i*nu*(z-z0)} + adf* u0*(z-z0)*e^{i*nu*(z-z0)}
!c  u1 is the solution to
!c  (\bTtil + nu*\bStil + nu^2*\bI).u1=i*adf*(\bStil + 2*nu*\bI).u0
!c  in practice, we absorb the adf factor into u1, then normalize
!c  (\bTtil + nu*\bStil + nu^2*\bI).(u1/adf)=i*(\bStil + 2*nu*\bI).u0
!c  since nu is the known eigenvalue of u0, the solution is easier
!c  form the matrices on either side
      znu=-eye*xnu(i,n)
      do ii=1,3
        do jj=1,3
          zq1(jj,ii)=dcmplx(ttl(jj,ii),0.d0)+znu*stl(jj,ii)
          zq2(jj,ii)=dcmplx(stl(jj,ii),0.d0)
        end do
        zq1(ii,ii)=zq1(ii,ii)+znu*znu
        zq2(ii,ii)=zq2(ii,ii)+2.d0*znu
      end do
!c  we wish to find the eigenvector of the near-defective matrix
!c  in the region where its eigenvectors are numerically unstable
!c   we explicitly calculate the eigenvector with smallest right-eigenvalue of
!c  (\bTtil + nu*\bStil + nu^2*\bI)=zq1
!c  copy into real, imag matrices
      do ii=1,3
        do jj=1,3
          q1r(jj,ii)=dreal(zq1(jj,ii))
          q1i(jj,ii)=dimag(zq1(jj,ii))
        end do
      end do      
!c  into eispack
      call cbal(3,3,q1r,q1i,low,igh,fv)
      call corth(3,3,low,igh,q1r,q1i,fv2,fv3)
      call comqr2(3,3,low,igh,fv2,fv3,q1r,q1i,wr,wi,q2r,q2i,ierr)
      if(ierr.ne.0) go to 400
      call cbabk2(3,3,low,igh,fv,3,q2r,q2i)
      amn=wr(1)**2+wi(1)**2
      ij=1
      do ii=2,3
        amm=wr(ii)**2+wi(ii)**2
        if(amm.lt.amn) then
          ij=ii
          amn=amm
        endif
      end do
      sum=0.d0
      do ii=1,3
        u0(ii)=dcmplx(q2r(ii,ij),q2i(ii,ij))
        sum=sum+zabs(u0(ii))**2
      end do
      sum=dsqrt(sum)
      do ii=1,3
        u0(ii)=u0(ii)/sum
      end do
!  assemble the ith stress-displacement vector
!  calculate the traction components, with i removed
      pp(1)=dcmplx(px,0.d0)
      pp(2)=z0
      pp(3)=znu
      pu=z0
      pw=z0
      uw=z0
      abcde=a-b+c-2.d0*d+2.d0*e
      bce=b-4.d0*c-4.d0*e
      de=d-e
      do ii=1,3
        pu=pu+pp(ii)*u0(ii)
        pw=pw+pp(ii)*w(ii,n)
        uw=uw+u0(ii)*w(ii,n)
      end do
      do ii=1,3
        ee(ii,i,n)= u0(ii)
        ee(ii+3,i,n)=w(ii,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c+2.d0*(pw*u0(3)+uw*pp(3))*e)
        ee(ii+3,i,n)=ee(ii+3,i,n)+pp(ii)*(u0(3)*de+2.d0*uw*w(3,n)*e)
        ee(ii+3,i,n)=ee(ii+3,i,n)+u0(ii)*(pp(3)*de+2.d0*pw*w(3,n)*e)
      end do
      ee(6,i,n)=ee(6,i,n)+pu*abcde+pw*uw*bce
!c  almost lastly, mult traction by i
      do ii=1,3
        ee(ii+3,i,n)=eye*ee(ii+3,i,n)
      end do
!c  extract u0 from ee(*,i,n) use it to calculate the additional traction terms
!c  and store in ee(*,j,n)
!c  additional traction terms involve gradient of (z-z0)
!c  so can be calculated from standard formulas with \bk=zhat
!c  we dont multiply by i
      pp(1)=z0
      pp(2)=z0
      pp(3)=z1
      pu=z0
      pw=z0
      uw=z0
      abcde=a-b+c-2.d0*d+2.d0*e
      bce=b-4.d0*c-4.d0*e
      de=d-e
      do ii=1,3
        u0(ii)=ee(ii,i,n)
        pu=pu+pp(ii)*u0(ii)
        pw=pw+pp(ii)*w(ii,n)
        uw=uw+u0(ii)*w(ii,n)
      end do
      do ii=1,3
        xee(ii)=w(ii,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c+2.d0*(pw*u0(3)+uw*pp(3))*e)
        xee(ii)=xee(ii)+pp(ii)*(u0(3)*de+2.d0*uw*w(3,n)*e)
        xee(ii)=xee(ii)+u0(ii)*(pp(3)*de+2.d0*pw*w(3,n)*e)
      end do
      xee(3)=xee(3)+pu*abcde+pw*uw*bce
!c  extract u0 from ee(*,i,n), mult by i*(\bStil + 2*nu*\bI), replace in u0
      do ii=1,3
        u0(ii)=z0
        do jj=1,3
          u0(ii)=u0(ii)+zq2(ii,jj)*ee(jj,i,n)
        end do
        u0(ii)=eye*u0(ii)
      end do
! 1002 format(3(2g14.6,3x))
!c  for znu NOT an eigenvalue,
!c  but rather the average of closely-space eigenvalues
!c  in this case, zq1 is nonsingular, and we just solve for u1 
      do ii=1,3
        yr(ii)=dreal(u0(ii))
        yi(ii)=dimag(u0(ii))
        do jj=1,3
          q1r(jj,ii)=dreal(zq1(jj,ii))
          q1i(jj,ii)=dimag(zq1(jj,ii))
        end do
      end do
      call csolve(3,q1r,q1i,xr,xi,yr,yi)
      do ii=1,3
        u1(ii)=dcmplx(xr(ii),xi(ii))
      end do
!c   End, different tactic
!c
!c  normalize
      sum=0.d0
      do ii=1,3
        sum=sum+zabs(u1(ii))**2
      end do
      sum=dsqrt(sum)
      do ii=1,3
        u1(ii)=u1(ii)/sum
      end do
!c  adf is the normalization constant
      adf=1.d0/sum
!c  calculate the traction
!c  and place the new stress-displacement vector in column j
!c  pp is the wavenumber vector, and first two components are already in place
      pp(1)=dcmplx(px,0.d0)
      pp(2)=z0
      pp(3)=znu
      pu=z0
      pw=z0
      uw=z0
      abcde=a-b+c-2.d0*d+2.d0*e
      bce=b-4.d0*c-4.d0*e
      de=d-e
      do ii=1,3
        pu=pu+pp(ii)*u1(ii)
        pw=pw+pp(ii)*w(ii,n)
        uw=uw+u1(ii)*w(ii,n)
      end do
      do ii=1,3
        ee(ii,j,n)=u1(ii)
        ee(ii+3,j,n)=w(ii,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c+2.d0*(pw*u1(3)+uw*pp(3))*e)
        ee(ii+3,j,n)=ee(ii+3,j,n)+pp(ii)*(u1(3)*de+2.d0*uw*w(3,n)*e)
        ee(ii+3,j,n)=ee(ii+3,j,n)+u1(ii)*(pp(3)*de+2.d0*pw*w(3,n)*e)
      end do
      ee(6,j,n)=ee(6,j,n)+pu*abcde+pw*uw*bce
!c  almost lastly, mult traction by i
!c  and add extra traction from (z-z0) term (not mult by i)
!c  TEST - mult xee by zero, see if it is important --- it IS important
      do ii=1,3
        ee(ii+3,j,n)=eye*ee(ii+3,j,n)+adf*xee(ii)
      end do
      return
 400  print *,'eispack error'
      stop
end

#endif

SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)

IMPLICIT NONE

INTEGER(kint) I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
REAL(kreal) A(6,6),SCALE(6)
REAL(kreal) C,F,G,R_balanc,S_balanc,B2,RADIX
REAL(kreal) DABS

LOGICAL NOCONV

!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
!     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).1

!     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
!     EIGENVALUES WHENEVER POSSIBLE.

!     ON INPUT:

!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;

!        N IS THE ORDER OF THE MATRIX;

!        A CONTAINS THE INPUT MATRIX TO BE BALANCED.

!     ON OUTPUT:

!        A CONTAINS THE BALANCED MATRIX;

!        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
!          IS EQUAL TO ZERO IF
!           (1) I IS GREATER THAN J AND
!           (2) J=1,...,LOW-1 OR I=IGH+1,...,N;

!        SCALE CONTAINS INFORMATION DETERMINING THE
!           PERMUTATIONS AND SCALING FACTORS USED.

!     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
!     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
!     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
!     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
!        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
!                 = D(J,J),      J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
!     THEN 1 TO LOW-1.

!     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.

!     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
!     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
!     K,L HAVE BEEN REVERSED.)

!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY

!     ------------------------------------------------------------------
!
!     :::::::::: RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
!                RADIX = 16.0D0 FOR LONG FORM ARITHMETIC
!                ON S360 ::::::::::
DATA RADIX/2/

B2 = RADIX * RADIX
K = 1
L = N

GO TO 100
!     :::::::::: IN-LINE PROCEDURE FOR ROW AND
!                COLUMN EXCHANGE ::::::::::
20 SCALE(M) = J
IF (J .EQ. M) GO TO 50
!
DO 30 I = 1, L
F = A(I,J)
A(I,J) = A(I,M)
A(I,M) = F
30 CONTINUE

DO 40 I = K, N
F = A(J,I)
A(J,I) = A(M,I)
A(M,I) = F
40 CONTINUE

50 GO TO (80,130), IEXC
!     :::::::::: SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                AND PUSH THEM DOWN ::::::::::
80 IF (L .EQ. 1) GO TO 280
L = L - 1
!     :::::::::: FOR J=L STEP -1 UNTIL 1 DO -- ::::::::::
100 DO 120 JJ = 1, L
J = L + 1 - JJ

DO 110 I = 1, L
IF (I .EQ. J) GO TO 110
IF (A(J,I) .NE. 0.0D0) GO TO 120
110    CONTINUE

M = L
IEXC = 1
GO TO 20
120 CONTINUE

GO TO 140
!     :::::::::: SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!                AND PUSH THEM LEFT ::::::::::
130 K = K + 1

140 DO 170 J = K, L

DO 150 I = K, L
IF (I .EQ. J) GO TO 150
IF (A(I,J) .NE. 0.0D0) GO TO 170
150    CONTINUE

M = K
IEXC = 2
GO TO 20
170 CONTINUE
!     :::::::::: NOW BALANCE THE SUBMATRIX IN ROWS K TO L ::::::::::
DO 180 I = K, L
180 SCALE(I) = 1.0D0
!     :::::::::: ITERATIVE LOOP FOR NORM REDUCTION ::::::::::

190 NOCONV = .FALSE.

DO 270 I = K, L
C = 0.0D0
R_balanc = 0.0D0

DO 200 J = K, L
!IF (J .EQ. I) GO TO 200
C = C + DABS(A(J,I))
R_balanc = R_balanc + DABS(A(I,J))
200    CONTINUE

!     :::::::::: GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ::::::::::
IF (C .EQ. 0.0D0 .OR. R_balanc .EQ. 0.0D0) GO TO 270
G = R_balanc / RADIX
F = 1.0D0
S_balanc = C + R_balanc

210    IF (C .GE. G) GO TO 220
F = F * RADIX
C = C * B2
GO TO 210
220    G = R_balanc * RADIX
230    IF (C .LT. G) GO TO 240
F = F / RADIX
C = C / B2
GO TO 230
!     :::::::::: NOW BALANCE ::::::::::
240    IF ((C + R_balanc) / F .GE. 0.95D0 * S_balanc) GO TO 270

G = 1.0D0 / F
SCALE(I) = SCALE(I) * F
NOCONV = .TRUE.

DO 250 J = K, N
250    A(I,J) = A(I,J) * G

DO 260 J = 1, L
260    A(J,I) = A(J,I) * F

270 CONTINUE

IF (NOCONV) GO TO 190

280 LOW = K
IGH = L

RETURN
!     :::::::::: LAST CARD OF BALANC ::::::::::
END



subroutine csolve(nn,a,ai,x_csolve6,xi_csolve6,y_csolve6,yi_csolve6, returnvalue)
!  solves the complex nxn system of equations a*x=y using gaussian elimination
!  and partial pivoting

implicit none

integer(kint), intent(in)           :: nn
real(kreal), intent(inout)          :: a(*), ai(*), x_csolve6(*), xi_csolve6(*), y_csolve6(*), yi_csolve6(*)

integer(kint)                       :: n, ip(1000), i, j

real(kreal), intent(inout) :: returnvalue

n=nn

if(n.gt.0) then
call clup(n,a,ai,ip)
else
n=-n
endif

call cbcktr(n,a,ai,x_csolve6,xi_csolve6,y_csolve6,yi_csolve6,ip, returnvalue)

return

end


subroutine clup(n,a,ai,ip)
!  finds lu decomp of a+i*ai using partial pivoting
!  pivoting sequence returned in ip(n)

implicit none

integer(kint), intent(in)          :: n
integer(kint), intent(inout)       :: ip(*)
real(kreal), intent(inout)         :: a(n,*), ai(n,*)

real(kreal)                        :: tol_clup, aam, aai, tem, b, bi
integer(kint)                      :: nm1, i, j, jm, i1, ipi, ipj, k

tol_clup=1.d-14

!  initialize permutation vector
do 50 i=1,n
50 ip(i)=i

nm1=n-1

if(n.eq.1) go to 700
do 100 i=1,nm1
aam=0.d0
do 200 j=i,n
aai=a(ip(j),i)**2+ai(ip(j),i)**2
if(aam.gt.aai) go to 200
aam=aai
ipi=ip(j)
jm=j
200 continue
if(aam.lt.tol_clup) go to 400
ip(jm)=ip(i)
ip(i)=ipi
i1=i+1
do 100 j=i1,n
ipj=ip(j)
! if victim index is already zero, dont bother to rub it out
tem=dabs(a(ipj,i))+dabs(ai(ipj,i))
if(tem.lt.tol_clup) go to 100
b=(a(ipj,i)*a(ipi,i)+ai(ipj,i)*ai(ipi,i))/aam
bi=(ai(ipj,i)*a(ipi,i)-a(ipj,i)*ai(ipi,i))/aam
a(ipj,i)=b
ai(ipj,i)=bi
do 500 k=i1,n
a(ipj,k)=a(ipj,k)-b*a(ipi,k)+bi*ai(ipi,k)
500 ai(ipj,k)=ai(ipj,k)-b*ai(ipi,k)-bi*a(ipi,k)
100 continue
700 continue
return
!c  400 print 101,aam,i
400 continue
!c  101 format('near-zero pivot ',e12.5,'  on column',i3)
!c      stop

end


subroutine cbcktr(n1,z_cbcktr6,zi_cbcktr6,dr,di,er,ei,ip, returnvalue)

implicit none

integer(kint), intent(in)           :: n1
real(kreal), intent(inout)          :: z_cbcktr6(n1,*), zi_cbcktr6(n1,*), dr(*), di(*), er(*), ei(*)
integer(kint), intent(inout)        :: ip(*)

integer(kint)                       :: i, j, i1, ii, ip1, iii, jjj
real(kreal)                         :: zkk1, zkk2, uii, dri, dii

real(kreal), intent(inout) :: returnvalue

!  back transform with unit lower triangular matrix
do 300 i=1,n1
dr(i)=er(ip(i))
300 di(i)=ei(ip(i))

if(n1.eq.1) go to 400
do 310 i=2,n1
i1=i-1
do 310 j=1,i1
zkk1=z_cbcktr6(ip(i),j)
zkk2=zi_cbcktr6(ip(i),j)
dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
310 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
400 continue
!  back transform with upper triangular matrix
do 320 ii=1,n1
i=n1+1-ii
ip1=i+1
uii=z_cbcktr6(ip(i),i)**2+zi_cbcktr6(ip(i),i)**2
if(i.eq.n1) go to 340
do 330 j=ip1,n1
zkk1=z_cbcktr6(ip(i),j)
zkk2=zi_cbcktr6(ip(i),j)
dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
330 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
340 dri=dr(i)
dii=di(i)
zkk1=z_cbcktr6(ip(i),i)
zkk2=zi_cbcktr6(ip(i),i)
dr(i)=(dri*zkk1+dii*zkk2)/uii
320 di(i)=(dii*zkk1-dri*zkk2)/uii

return
end


subroutine solve(nn,a_solve,x_solve,y_solve, returnvalue)

!  solves the nxn system of equations a*x=y using gaussian elimination
!  and partial pivoting
!  if n<0 the lu decomposition is already done
!  note that the matrix a is modified

implicit none

real(kreal), intent(inout)      :: a_solve(*)

real(kreal), intent(inout)      :: x_solve(*), y_solve(*)
integer(kint), intent(in)       :: nn
integer(kint)                   :: ip(1000)

real(kreal), intent(inout)      :: returnvalue

integer(kint)                   :: n, i, j, iii, jjj

n=nn

if(n.gt.0)then

    call lup(n,a_solve,ip)

else

    n=-n

endif

!      print 101,((i,j,a(i,j),j=1,n),i=1,n)
!  101 format(' a(',2i2,')=',e15.5)
!      type 102,(ip(i),i=1,n)
!  102 format(10i5)

call bcktr(n,a_solve,x_solve,y_solve,ip)

return
end


subroutine lup(n,a_lup,ip)

implicit none

!  finds lu decomp of a using partial pivoting
!  output in c (upper triangle/unit lower triangle) and
!  pivoting sequence returned in ip(n)

integer(kint), intent(in)     :: n
integer(kint), intent(inout)  :: ip(*)
real(kreal), intent(inout)    :: a_lup(n,*)

real(kreal)                   :: tol_lup, aam, aai, tem, b
integer(kint)                 :: nm1, i, j, jm, ipi, i1, ipj,k

tol_lup=1.d-14
do 50 i=1,n
50 ip(i)=i
nm1=n-1
if(n.eq.1) go to 700
do 100 i=1,nm1
aam=0.d0
do 200 j=i,n
aai=a_lup(ip(j),i)**2
if(aam.gt.aai) go to 200
aam=aai
ipi=ip(j)
jm=j
200 continue
if(aam.lt.tol_lup) go to 400
ip(jm)=ip(i)
ip(i)=ipi
i1=i+1
do 100 j=i1,n
ipj=ip(j)
!  if victim index is already zero, dont bother to rub it out
tem=dabs(a_lup(ipj,i))
if(tem.lt.tol_lup) go to 100
b=(a_lup(ipj,i)*a_lup(ipi,i))/aam
a_lup(ipj,i)=b
do 500 k=i1,n
a_lup(ipj,k)=a_lup(ipj,k)-b*a_lup(ipi,k)
500 continue
100 continue
700 continue
return
!  400 print 101,aam,i
400 continue
!  101 format('near-zero pivot ',e12.5,'  on column',i3)
!      stop


end

subroutine bcktr(n1,z_bcktr,dr,er,ip)

implicit none

!  performs backtransform on input vector er - 'y'
!  to find solution dr - 'x'

real(kreal), intent(in)      :: z_bcktr(n1, *)
real(kreal), intent(inout)   :: dr(*), er(*)
integer(kint), intent(inout) :: ip(*), n1

integer(kint) i, i1, ii, ip1, j

!  back transform with unit lower triangular matrix

do 300 i=1,n1
300 dr(i)=er(ip(i))
if(n1.eq.1) go to 400
do 310 i=2,n1
i1=i-1
do 310 j=1,i1
310 dr(i)=dr(i)-z_bcktr(ip(i),j)*dr(j)
400 continue

!  back transform with upper triangular matrix
do 320 ii=1,n1
i=n1+1-ii
ip1=i+1
if(i.eq.n1) go to 320
do 330 j=ip1,n1
330 dr(i)=dr(i)-z_bcktr(ip(i),j)*dr(j)
320 dr(i)=dr(i)/z_bcktr(ip(i),i)

return

end

SUBROUTINE elmhes(NM,N,LOW,IGH,A,INT)

implicit none

Integer(kint)                           :: I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
REAL(kreal), intent(inout)              :: A(NM,N)
REAL(kreal)                             :: X_elmhes,Y_elmhes

INTEGER(kint), intent(inout)            :: INT(IGH)


!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMHES,
!C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!C
!C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
!C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
!C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
!C     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.
!C
!C     ON INPUT:
!C
!C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!C          DIMENSION STATEMENT;
!C
!C        N IS THE ORDER OF THE MATRIX;
!C
!C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
!C          SET LOW=1, IGH=N;
!C
!C        A CONTAINS THE INPUT MATRIX.
!C
!C     ON OUTPUT:
!C
!C        A CONTAINS THE HESSENBERG MATRIX.  THE MULTIPLIERS
!C          WHICH WERE USED IN THE REDUCTION ARE STORED IN THE
!C          REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX;
!C
!C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
!C          INTERCHANGED IN THE REDUCTION.
!C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
!C
!C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------

LA = IGH - 1
KP1 = LOW + 1
IF (LA .LT. KP1) GO TO 200

DO 180 M = KP1, LA
MM1 = M - 1
X_elmhes = 0.0D0
I = M

DO 100 J = M, IGH
IF (DABS(A(J,MM1)) .LE. DABS(X_elmhes)) GO TO 100
X_elmhes = A(J,MM1)
I = J
100    CONTINUE

INT(M) = I
IF (I .EQ. M) GO TO 130
!    :::::::::: INTERCHANGE ROWS AND COLUMNS OF A ::::::::::
DO 110 J = MM1, N
Y_elmhes = A(I,J)
A(I,J) = A(M,J)
A(M,J) = Y_elmhes
110    CONTINUE

DO 120 J = 1, IGH
Y_elmhes = A(J,I)
A(J,I) = A(J,M)
A(J,M) = Y_elmhes
120    CONTINUE
!    :::::::::: END INTERCHANGE ::::::::::
130    IF (X_elmhes .EQ. 0.0D0) GO TO 180
MP1 = M + 1

DO 160 I = MP1, IGH
Y_elmhes = A(I,MM1)
IF (Y_elmhes .EQ. 0.0D0) GO TO 160
Y_elmhes = Y_elmhes / X_elmhes
A(I,MM1) = Y_elmhes

DO 140 J = M, N
140       A(I,J) = A(I,J) - Y_elmhes * A(M,J)

DO 150 J = 1, IGH
150       A(J,M) = A(J,M) + Y_elmhes * A(J,I)

160    CONTINUE

180 CONTINUE

200 RETURN
!    :::::::::: LAST CARD OF ELMHES ::::::::::
END

SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z_elt)

implicit none

INTEGER(kint)                 :: I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
REAL(kreal), intent(in)       :: A(6,6)
REAL(kreal), intent(inout)    :: Z_elt(6,6)
INTEGER(kint)                 :: INT(6)

!C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMTRANS,
!C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!C
!C     THIS SUBROUTINE ACCUMULATES THE STABILIZED ELEMENTARY
!C     SIMILARITY TRANSFORMATIONS USED IN THE REDUCTION OF A
!C     REAL GENERAL MATRIX TO UPPER HESSENBERG FORM BY  ELMHES.
!C
!C     ON INPUT:
!C
!C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!C          DIMENSION STATEMENT;
!C
!C        N IS THE ORDER OF THE MATRIX;
!C
!C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
!C          SET LOW=1, IGH=N;
!C
!C        A CONTAINS THE MULTIPLIERS WHICH WERE USED IN THE
!C          REDUCTION BY  ELMHES  IN ITS LOWER TRIANGLE
!C          BELOW THE SUBDIAGONAL;
!C
!C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
!C          INTERCHANGED IN THE REDUCTION BY  ELMHES.
!C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
!C
!C     ON OUTPUT:
!C
!C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!C          REDUCTION BY  ELMHES.
!C
!C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C
!C     ------------------------------------------------------------------
!C
!     :::::::::: INITIALIZE Z TO IDENTITY MATRIX ::::::::::
DO 80 I = 1, N

DO 60 J = 1, N
60    Z_elt(I,J) = 0.0D0

Z_elt(I,I) = 1.0D0
80 CONTINUE

KL = IGH - LOW - 1
IF (KL .LT. 1) GO TO 200
!C     :::::::::: FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ::::::::::
DO 140 MM = 1, KL
MP = IGH - MM
MP1 = MP + 1

DO 100 I = MP1, IGH
100    Z_elt(I,MP) = A(I,MP-1)

I = INT(MP)
IF (I .EQ. MP) GO TO 140

DO 130 J = MP, IGH
Z_elt(MP,J) = Z_elt(I,J)
Z_elt(I,J) = 0.0D0
130    CONTINUE

Z_elt(I,MP) = 1.0D0
140 CONTINUE

200 RETURN
!C     :::::::::: LAST CARD OF ELTRAN ::::::::::
END

SUBROUTINE hqr2(NM,N,LOW,IGH,H,WR_hqr2,WI_hqr2,Z_hqr2,IERR)

implicit none

INTEGER(kint)          :: I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,IGH,ITS,LOW,MP2,ENM2,IERR
REAL(kreal)            :: H(6,6),WR_hqr2(6),WI_hqr2(6),Z_hqr2(6,6)
REAL(kreal)            :: P,Q,R_hqr2,S_hqr2,T_hqr2,W_hqr2,X_hqr2,Y_hqr2,RA,SA,VI,VR,ZZ,NORM,MACHEP
REAL(kreal)            :: DSQRT,DABS,DSIGN
INTEGER(kint)          :: MIN0

LOGICAL                :: NOTLAS

COMPLEX(kreal)         :: Z3
COMPLEX(kreal)         :: DCMPLX
REAL(kreal)            :: DREAL,DIMAG

!C     :::::::::: STATEMENT FUNCTIONS ENABLE EXTRACTION OF REAL AND
!C     IMAGINARY PARTS OF DOUBLE PRECISION COMPLEX NUMBERS ::::::::::

DREAL(Z3) = Z3
DIMAG(Z3) = (0.0D0,-1.0D0) * Z3

!C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
!C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!C
!C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
!C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
!C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
!C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
!C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
!C
!C     ON INPUT:
!C
!C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!C          DIMENSION STATEMENT;
!C
!C        N IS THE ORDER OF THE MATRIX;
!C
!C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
!C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
!C          SET LOW=1, IGH=N;
!C
!C        H CONTAINS THE UPPER HESSENBERG MATRIX;
!C
!C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
!C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
!C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
!C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
!C          IDENTITY MATRIX.
!C
!C     ON OUTPUT:
!C
!C        H HAS BEEN DESTROYED;
!C
!C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
!C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
!C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
!C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
!C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
!C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
!C          FOR INDICES IERR+1,...,N;
!C
!C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
!C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
!C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
!C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
!C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
!C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
!C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND;
!C
!C        IERR IS SET TO
!C          ZERO       FOR NORMAL RETURN,
!C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!C                     DETERMINED AFTER 30 ITERATIONS.
!C
!C     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
!C     PROCEDURE CDIV BY COMPLEX DIVISION.
!C
!C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C
!C     ------------------------------------------------------------------
!C
!C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
!C                ON S360 ::::::::::
!C  FOR F_FLOATING DEC FORTRAN
!C      DATA MACHEP/1.1D-16/
!C  FOR G_FLOATING DEC FORTRAN

DATA MACHEP/1.25D-15/

IERR = 0
NORM = 0.0D0
K = 1
!C     :::::::::: STORE ROOTS ISOLATED BY BALANC
!C                AND COMPUTE MATRIX NORM ::::::::::
DO 50 I = 1, N

DO 40 J = K, N
40    NORM = NORM + DABS(H(I,J))

K = I
IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
WR_hqr2(I) = H(I,I)
WI_hqr2(I) = 0.0D0
50 CONTINUE

EN = IGH
T = 0.0D0
!C     :::::::::: SEARCH FOR NEXT EIGENVALUES ::::::::::
60 IF (EN .LT. LOW) GO TO 340
ITS = 0
NA = EN - 1
ENM2 = NA - 1
!C     :::::::::: LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!C                FOR L=EN STEP -1 UNTIL LOW DO -- ::::::::::
70 DO 80 LL = LOW, EN
L = EN + LOW - LL
IF (L .EQ. LOW) GO TO 100
S_hqr2 = DABS(H(L-1,L-1)) + DABS(H(L,L))
IF (S_hqr2 .EQ. 0.0D0) S_hqr2 = NORM
IF (DABS(H(L,L-1)) .LE. MACHEP * S_hqr2) GO TO 100
80 CONTINUE
!C     :::::::::: FORM SHIFT ::::::::::
100 X_hqr2 = H(EN,EN)
IF (L .EQ. EN) GO TO 270
Y_hqr2 = H(NA,NA)
W_hqr2 = H(EN,NA) * H(NA,EN)
IF (L .EQ. NA) GO TO 280
IF (ITS .EQ. 30) GO TO 1000
IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
!C     :::::::::: FORM EXCEPTIONAL SHIFT ::::::::::
T_hqr2 = T_hqr2 + X_hqr2

DO 120 I = LOW, EN
120 H(I,I) = H(I,I) - X_hqr2

S_hqr2 = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
X_hqr2= 0.75D0 * S_hqr2
Y_hqr2 = X_hqr2
W_hqr2 = -0.4375D0 * S_hqr2 * S_hqr2
130 ITS = ITS + 1
!C     :::::::::: LOOK FOR TWO CONSECUTIVE SMALL
!C                SUB-DIAGONAL ELEMENTS.
!C                FOR M=EN-2 STEP -1 UNTIL L DO -- ::::::::::
DO 140 MM = L, ENM2
M = ENM2 + L - MM
ZZ = H(M,M)
R_hqr2 = X_hqr2 - ZZ
S_hqr2 = Y_hqr2 - ZZ
P = (R_hqr2 * S_hqr2 - W_hqr2) / H(M+1,M) + H(M,M+1)
Q = H(M+1,M+1) - ZZ - R_hqr2 - S_hqr2
R_hqr2 = H(M+2,M+1)
S_hqr2 = DABS(P) + DABS(Q) + DABS(R_hqr2)
P = P / S_hqr2
Q = Q / S_hqr2
R_hqr2 = R_hqr2 / S_hqr2
IF (M .EQ. L) GO TO 150
IF (DABS(H(M,M-1)) * (DABS(Q) + DABS(R_hqr2)) .LE. MACHEP * DABS(P) * (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150
140 CONTINUE

150 MP2 = M + 2

DO 160 I = MP2, EN
H(I,I-2) = 0.0D0
IF (I .EQ. MP2) GO TO 160
H(I,I-3) = 0.0D0
160 CONTINUE
!C     :::::::::: DOUBLE QR STEP INVOLVING ROWS L TO EN AND
!C                COLUMNS M TO EN ::::::::::
DO 260 K = M, NA
NOTLAS = K .NE. NA
IF (K .EQ. M) GO TO 170
P = H(K,K-1)
Q = H(K+1,K-1)
R_hqr2 = 0.0D0
IF (NOTLAS) R_hqr2 = H(K+2,K-1)
X_hqr2 = DABS(P) + DABS(Q) + DABS(R_hqr2)
IF (X_hqr2 .EQ. 0.0D0) GO TO 260
P = P / X_hqr2
Q = Q / X_hqr2
R_hqr2 = R_hqr2 / X_hqr2
170    S_hqr2 = DSIGN(DSQRT(P*P+Q*Q+R_hqr2*R_hqr2),P)
IF (K .EQ. M) GO TO 180
H(K,K-1) = -S_hqr2 * X_hqr2
GO TO 190
180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
190    P = P + S_hqr2
X_hqr2 = P / S_hqr2
Y_hqr2 = Q / S_hqr2
ZZ = R_hqr2 / S_hqr2
Q = Q / P
R_hqr2 = R_hqr2 / P
!C     :::::::::: ROW MODIFICATION ::::::::::
DO 210 J = K, N
P = H(K,J) + Q * H(K+1,J)
IF (.NOT. NOTLAS) GO TO 200
P = P + R_hqr2 * H(K+2,J)
H(K+2,J) = H(K+2,J) - P * ZZ
200       H(K+1,J) = H(K+1,J) - P * Y_hqr2
H(K,J) = H(K,J) - P * X_hqr2
210    CONTINUE

J = MIN0(EN,K+3)
!C     :::::::::: COLUMN MODIFICATION ::::::::::
DO 230 I = 1, J
P = X_hqr2 * H(I,K) + Y_hqr2 * H(I,K+1)
IF (.NOT. NOTLAS) GO TO 220
P = P + ZZ * H(I,K+2)
H(I,K+2) = H(I,K+2) - P * R_hqr2
220       H(I,K+1) = H(I,K+1) - P * Q
H(I,K) = H(I,K) - P
230    CONTINUE
!C     :::::::::: ACCUMULATE TRANSFORMATIONS ::::::::::
DO 250 I = LOW, IGH
P = X_hqr2 * Z_hqr2(I,K) + Y_hqr2 * Z_hqr2(I,K+1)
IF (.NOT. NOTLAS) GO TO 240
P = P + ZZ * Z_hqr2(I,K+2)
Z_hqr2(I,K+2) = Z_hqr2(I,K+2) - P * R_hqr2
240       Z_hqr2(I,K+1) = Z_hqr2(I,K+1) - P * Q
Z_hqr2(I,K) = Z_hqr2(I,K) - P
250    CONTINUE

260 CONTINUE

GO TO 70
!C     :::::::::: ONE ROOT FOUND ::::::::::
270 H(EN,EN) = X_hqr2 + T_hqr2
WR_hqr2(EN) = H(EN,EN)
WI_hqr2(EN) = 0.0D0
EN = NA
GO TO 60
!C     :::::::::: TWO ROOTS FOUND ::::::::::
280 P = (Y_hqr2 - X_hqr2) / 2.0D0
Q = P * P + W_hqr2
ZZ = DSQRT(DABS(Q))
H(EN,EN) = X_hqr2 + T_hqr2
X_hqr2 = H(EN,EN)
H(NA,NA) = Y_hqr2 + T_hqr2
IF (Q .LT. 0.0D0) GO TO 320
!C     :::::::::: REAL PAIR ::::::::::
ZZ = P + DSIGN(ZZ,P)
WR_hqr2(NA) = X_hqr2 + ZZ
WR_hqr2(EN) = WR_hqr2(NA)
IF (ZZ .NE. 0.0D0) WR_hqr2(EN) = X_hqr2 - W_hqr2 / ZZ
WI_hqr2(NA) = 0.0D0
WI_hqr2(EN) = 0.0D0
X_hqr2 = H(EN,NA)
S_hqr2 = DABS(X_hqr2) + DABS(ZZ)
P = X_hqr2 / S_hqr2
Q = ZZ / S_hqr2
R_hqr2 = DSQRT(P*P+Q*Q)
P = P / R_hqr2
Q = Q / R_hqr2
!C     :::::::::: ROW MODIFICATION ::::::::::
DO 290 J = NA, N
ZZ = H(NA,J)
H(NA,J) = Q * ZZ + P * H(EN,J)
H(EN,J) = Q * H(EN,J) - P * ZZ
290 CONTINUE
!C     :::::::::: COLUMN MODIFICATION ::::::::::
DO 300 I = 1, EN
ZZ = H(I,NA)
H(I,NA) = Q * ZZ + P * H(I,EN)
H(I,EN) = Q * H(I,EN) - P * ZZ
300 CONTINUE
!C     :::::::::: ACCUMULATE TRANSFORMATIONS ::::::::::
DO 310 I = LOW, IGH
ZZ = Z_hqr2(I,NA)
Z_hqr2(I,NA) = Q * ZZ + P * Z_hqr2(I,EN)
Z_hqr2(I,EN) = Q * Z_hqr2(I,EN) - P * ZZ
310 CONTINUE

GO TO 330
!C     :::::::::: COMPLEX PAIR ::::::::::
320 WR_hqr2(NA) = X_hqr2 + P
WR_hqr2(EN) = X_hqr2 + P
WI_hqr2(NA) = ZZ
WI_hqr2(EN) = -ZZ
330 EN = ENM2
GO TO 60
!C     :::::::::: ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
!C                VECTORS OF UPPER TRIANGULAR FORM ::::::::::
340 IF (NORM .EQ. 0.0D0) GO TO 1001
!C     :::::::::: FOR EN=N STEP -1 UNTIL 1 DO -- ::::::::::
DO 800 NN = 1, N
EN = N + 1 - NN
P = WR_hqr2(EN)
Q = WI_hqr2(EN)
NA = EN - 1
IF (Q) 710, 600, 800
!C     :::::::::: REAL VECTOR ::::::::::
600    M = EN
H(EN,EN) = 1.0D0
IF (NA .EQ. 0) GO TO 800
!C     :::::::::: FOR I=EN-1 STEP -1 UNTIL 1 DO -- ::::::::::
DO 700 II = 1, NA
I = EN - II
W_hqr2 = H(I,I) - P
R_hqr2 = H(I,EN)
IF (M .GT. NA) GO TO 620

DO 610 J = M, NA
610       R_hqr2 = R_hqr2 + H(I,J) * H(J,EN)

620       IF (WI_hqr2(I) .GE. 0.0D0) GO TO 630
ZZ = W_hqr2
S_hqr2 = R_hqr2
GO TO 700
630       M = I
IF (WI_hqr2(I) .NE. 0.0D0) GO TO 640
T_hqr2 = W_hqr2
IF (W_hqr2 .EQ. 0.0D0) T_hqr2 = MACHEP * NORM
H(I,EN) = -R_hqr2 / T_hqr2
GO TO 700
!C     :::::::::: SOLVE REAL EQUATIONS ::::::::::
640       X_hqr2 = H(I,I+1)
Y_hqr2 = H(I+1,I)
Q = (WR_hqr2(I) - P) * (WR_hqr2(I) - P) + WI_hqr2(I) * WI_hqr2(I)
T_hqr2 = (X_hqr2 * S_hqr2 - ZZ * R_hqr2) / Q
H(I,EN) = T_hqr2
IF (DABS(X_hqr2) .LE. DABS(ZZ)) GO TO 650
H(I+1,EN) = (-R_hqr2 - W_hqr2 * T_hqr2) / X_hqr2
GO TO 700
650       H(I+1,EN) = (-S_hqr2 - Y_hqr2 * T_hqr2) / ZZ
700    CONTINUE
!C     :::::::::: END REAL VECTOR ::::::::::
GO TO 800
!C     :::::::::: COMPLEX VECTOR ::::::::::
710    M = NA
!C     :::::::::: LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
!C                EIGENVECTOR MATRIX IS TRIANGULAR ::::::::::
IF (DABS(H(EN,NA)) .LE. DABS(H(NA,EN))) GO TO 720
H(NA,NA) = Q / H(EN,NA)
H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
GO TO 730
720    Z3 = DCMPLX(0.0D0,-H(NA,EN)) / DCMPLX(H(NA,NA)-P,Q)
H(NA,NA) = DREAL(Z3)
H(NA,EN) = DIMAG(Z3)
730    H(EN,NA) = 0.0D0
H(EN,EN) = 1.0D0
ENM2 = NA - 1
IF (ENM2 .EQ. 0) GO TO 800
!C     :::::::::: FOR I=EN-2 STEP -1 UNTIL 1 DO -- ::::::::::
DO 790 II = 1, ENM2
I = NA - II
W_hqr2 = H(I,I) - P
RA = 0.0D0
SA = H(I,EN)

DO 760 J = M, NA
RA = RA + H(I,J) * H(J,NA)
SA = SA + H(I,J) * H(J,EN)
760       CONTINUE

IF (WI(I) .GE. 0.0D0) GO TO 770
ZZ = W_hqr2
R_hqr2 = RA
S_hqr2 = SA
GO TO 790
770       M = I
IF (WI_hqr2(I) .NE. 0.0D0) GO TO 780
Z3 = DCMPLX(-RA,-SA) / DCMPLX(W_hqr2,Q)
H(I,NA) = DREAL(Z3)
H(I,EN) = DIMAG(Z3)
GO TO 790
!C     :::::::::: SOLVE COMPLEX EQUATIONS ::::::::::
780       X_hqr2 = H(I,I+1)
Y_hqr2 = H(I+1,I)
VR = (WR_hqr2(I) - P) * (WR_hqr2(I) - P) + WI_hqr2(I) * WI_hqr2(I) - Q * Q
VI = (WR_hqr2(I) - P) * 2.0D0 * Q
IF (VR .EQ. 0.0D0 .AND. VI .EQ. 0.0D0) VR = MACHEP * NORM * (DABS(W_hqr2) + DABS(Q) + DABS(X_hqr2) + DABS(Y_hqr2) + DABS(ZZ))
Z3 = DCMPLX(X_hqr2*R_hqr2-ZZ*RA+Q*SA,X_hqr2*S_hqr2-ZZ*SA-Q*RA) / DCMPLX(VR,VI)
H(I,NA) = DREAL(Z3)
H(I,EN) = DIMAG(Z3)
IF (DABS(X_hqr2) .LE. DABS(ZZ) + DABS(Q)) GO TO 785
H(I+1,NA) = (-RA - W_hqr2 * H(I,NA) + Q * H(I,EN)) / X_hqr2
H(I+1,EN) = (-SA - W_hqr2 * H(I,EN) - Q * H(I,NA)) / X_hqr2
GO TO 790
785       Z3 = DCMPLX(-R_hqr2-Y_hqr2*H(I,NA),-S_hqr2-Y_hqr2*H(I,EN)) / DCMPLX(ZZ,Q)
H(I+1,NA) = DREAL(Z3)
H(I+1,EN) = DIMAG(Z3)
790    CONTINUE
!C     :::::::::: END COMPLEX VECTOR ::::::::::
800 CONTINUE
!     :::::::::: END BACK SUBSTITUTION.
!                VECTORS OF ISOLATED ROOTS ::::::::::
DO 840 I = 1, N
IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840

DO 820 J = I, N
820    Z_hqr2(I,J) = H(I,J)

840 CONTINUE
!C     :::::::::: MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!C                VECTORS OF ORIGINAL FULL MATRIX.
!C                FOR J=N STEP -1 UNTIL LOW DO -- ::::::::::
DO 880 JJ = LOW, N
J = N + LOW - JJ
M = MIN0(J,IGH)

DO 880 I = LOW, IGH
ZZ = 0.0D0

DO 860 K = LOW, M
860       ZZ = ZZ + Z_hqr2(I,K) * H(K,J)

Z_hqr2(I,J) = ZZ
880 CONTINUE

GO TO 1001
!C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
!C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
1000 IERR = EN
1001 RETURN
!C     :::::::::: LAST CARD OF HQR2 ::::::::::


END


SUBROUTINE balbak(NM,N,LOW,IGH,SCALE,M,Z_balbak)

implicit none

INTEGER(kint)                          :: I,J,K,M,N,II,NM,IGH,LOW
REAL(kreal)                            :: SCALE(6)
real(kreal), intent(inout)             :: Z_balbak(6,6)
REAL(kreal)                            :: S_balbak

!C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALBAK,
!C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
!C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
!C
!C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
!C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!C     BALANCED MATRIX DETERMINED BY  BALANC.
!C
!C     ON INPUT:
!C
!C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!C          DIMENSION STATEMENT;
!C
!C        N IS THE ORDER OF THE MATRIX;
!C
!C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC;
!C
!C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
!C          AND SCALING FACTORS USED BY  BALANC;
!C
!C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED;
!C
!C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
!C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
!C
!C     ON OUTPUT:
!C
!C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
!C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
!C
!C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!C
!C     ------------------------------------------------------------------

IF (M .EQ. 0) GO TO 200
IF (IGH .EQ. LOW) GO TO 120

DO 110 I = LOW, IGH
S_balbak = SCALE(I)
!     :::::::::: LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!                IF THE FOREGOING STATEMENT IS REPLACED BY
!                S=1.0D0/SCALE(I). ::::::::::
DO 100 J = 1, M
100    Z_balbak(I,J) = Z_balbak(I,J) * S_balbak

110 CONTINUE
!     ::::::::: FOR I=LOW-1 STEP -1 UNTIL 1,
!               IGH+1 STEP 1 UNTIL N DO -- ::::::::::
120 DO 140 II = 1, N
I = II
IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
IF (I .LT. LOW) I = LOW - II
K = SCALE(I)
IF (K .EQ. I) GO TO 140

DO 130 J = 1, M
S_balbak = Z_balbak(I,J)
Z_balbak(I,J) = Z_balbak(K,J)
Z_balbak(K,J) = S_balbak
130    CONTINUE

140 CONTINUE

200 RETURN
!     :::::::::: LAST CARD OF BALBAK ::::::::::
END



END MODULE anirec_module

