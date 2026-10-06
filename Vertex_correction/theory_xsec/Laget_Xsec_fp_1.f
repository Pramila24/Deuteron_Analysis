c------------------------------------------------------------------------------


c
c calculate interpolated cross sections for a given
c kinematics using j.m.laget's reponse function
c 
c this code has been extraced from mceep
c same as Laget_Xsec.f but one can select the calculation of FSi
c with the parametere do_fsi: = 1 yes, =0 no
c------------------------------------------------------------------------------
c
c
      subroutine init_laget( datadir, do_fsi, interpol_type, 
     >     save_data, use_binary_file)

      integer n_dat_dir, do_fsi, interpol_type
      integer save_data, use_binary_file
      character*100 dat_dir
      character*100 datadir

      logical save_grid, use_binary
      common/fc__control_c/ save_grid, use_binary

      common/fc__datdir_c/ dat_dir
      common/fc__datdir_i/ n_dat_dir

      real*8 xm_e,xm_p,xm_n,xm_d,xmd,spn,dpn,xme2,xmn2,xmp2,xmd2,spn2
      common/fc__some_mas / dpn,xmd,spn2,xmd2,xmn2

c laget data grid

      common/fc__grid_sig / sig_l,sig_t,sig_lt,sig_tt
      common/fc__laget_mod/ laget_intp,laget_pwia,laget_fsi,laget_mec

      real*8 sig_l(100,500,91),sig_t(100,500,91)
      real*8 sig_lt(100,500,91), sig_tt(100,500,91)


      data  xm_e / 510.99890d-03 /,  xm_p / 938.27200d+00 /, 
     +      xm_n / 939.56533d+00 /,  xm_d / 187.56128d+01 /


c
c laget calculation selectors
c
      integer laget_intp,laget_pwia,laget_fsi,laget_mec
      integer error
      integer i1, i2
c
c
c ------------------------------------------------------------------------
c       d(e,e'p)n laget unpolarized response functions
c ------------------------------------------------------------------------
c
      integer lin, log, scat_neutron, scat_proton, scat_both
      integer fsi, no_fsi, mec, no_mec

      parameter (lin=1)
      parameter (log = 2)
      parameter (scat_neutron = 0)
      parameter (scat_proton = 1)
      parameter (scat_both = 2)
      parameter (no_fsi = 0)
      parameter (fsi = 1)
      parameter (no_mec = 0)
      parameter (mec = 1)

*---- particle masses and relevant combinations                             

      xmd = xm_d
*-    dimension 1 combinations
      spn = xm_n + xm_p
      dpn = xm_n - xm_p
*-    dimension 2 combinations                                     
      xme2 = xm_e * xm_e                                            
      xmd2 = xm_d * xm_d                                                       
      xmp2 = xm_p * xm_p
      xmn2 = xm_n * xm_n
      spn2 =  spn * spn                                                      

      if (save_data .eq. 1) then
         save_grid = .True.
      else
         save_grid = .False.
      endif

      if (use_binary_file .eq. 1) then
         use_binary = .True.
      else
         use_binary = .False.
      endif

c make sure there is non nonsense
      if (save_grid) use_binary = .False.
      if (use_binary) save_grid = .False.

c
c standard configuration
c
c use log. interpolation (better)
c
      laget_intp = interpol_type
      laget_pwia = scat_both
      laget_fsi = do_fsi
      laget_mec = no_mec
c
c setup data directory name
c      
      call fs__nospac(datadir, i1, i2)
      dat_dir = datadir(i1:i2)
      n_dat_dir = len(dat_dir)

      write(*,'(a)') ' start loading the data grid...'  

      

      call fs__grid_load(sig_l,sig_t,sig_lt,sig_tt, error)            
      if (error .ne. 0) then
         write(*,'(a)') ' problem reading of data grid load !'
         write(*,'(a)') ' I would not continue !'
         return
      endif
      write(*,'(a)') ' end of data grid load !'  
      
      
c     
      return
      end




c------------------------------------------------------------------------------
c 
c this file is also callable from python
c


      real*8 function get_sigma_laget(e0_i,qmu2,omega,theta_cm,phi_x)
C------------------------------------------------------------------------------
c       purpose:
c               get (e,e'n) cross sections and polarizations
c               according to choices made in subroutine phys_choice.
c
c               coincidence cross sections (sigma_eep) are returned
c               with the following units:
c                      fm^2 sr^-2 mev^-1            (bound state)
c                      fm^2 sr^-2 mev^-1 (mev/c)^-1 (continuum)
c               note that the continuum cross section should be
c               differential in the hadron momentum (not kinetic energy)
c               since it is the momentum which is randomly sampled.
c               the recoil factor for the continuum case is equal to
c               the ejectile total energy divided by the momentum.
c               division by this factor converts the cross section
c               from being differential in energy (as, for example
c               is the case for the deforest cc1 cross section) to
c               differential in momentum.
c ---------------------------------------------------------------------
c
c
      implicit none
c
c arguments
c
      real*8 e0_i,qmu2,omega,theta_cm,phi_x
c     
c return value
c
      real*8 sigma_eep

c grid
      common/fc__grid_sig / sig_l,sig_t,sig_lt,sig_tt
c
c laget grid
      real*8 sig_l(100,500,91),sig_t(100,500,91)
      real*8 sig_lt(100,500,91), sig_tt(100,500,91)

c
      real*8 hbarc, pi
      parameter (hbarc = 197.3286d0)
      parameter (pi=3.14159265359d0) 
c
c ------------------------------------------------------------------------
c       laget unpolarized response functions (j.-m. laget)
c ------------------------------------------------------------------------
c
      call fs__laget_xsec(sigma_eep,e0_i,qmu2,omega,theta_cm,phi_x,
     +     sig_l,sig_t,sig_lt,sig_tt)
      get_sigma_laget = sigma_eep * 1.d-4  ! ub -> fm^2
      return 

      end
c
c------------------------------------------------------------------------------
c    

      subroutine fs__laget_xsec(xsec,p_bea,g_mas,g_ene,p_tgs,p_phi,
     +                      sig_l,sig_t,sig_lt,sig_tt)

c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
c
c  author: e. voutier
c  date: october 2005
c  purpose: determine d(e,e'p) cross section at current kinematics 
c  
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c  
c  calculate the d(e,e'p) cross section for a selected kinematics from an 
c  input data  grid via
c  the interpolation  procedure specified by the user via the laget_intp  
c  parameter (1 for linear, 2 for logarythmic in recoil momentum).  
c 
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  input variables : 
c
c        p_bea       3-momentum of the electron beam (mev/c)         
c        g_mas       quadrimomentum transfer (mev^2)
c        g_ene       energy transfer (mev)
c        p_tgs       proton angle in the center of mass frame (rd)
c        p_phi       out-of-plane angle of the proton (rd)
c        sig_l       array of the longitudinal cross sections
c        sig_t       array of the transverse cross sections
c        sig_lt      array of the longitudinal-transverse cross sections
c        sig_tt      array of the transverse-transverse cross sections
c
c  passed via common :
c
c        laget_intp  interpolation index
c
c  output variables :
c 
c        xsec    actual d(e,e'p) cross section (ub.mev-1.sr-2)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

      implicit none
      
      include 'grid_par.inc'

      integer   laget_ps_fail,laget_grid_fail
      integer*4 laget_intp,laget_pwia,laget_fsi,laget_mec

      real*8 ff__intrpol,ff__intrlog    
            
      real*8 xsec,p_bea,g_mas,g_ene,p_tgs,p_phi
       
      real*8  sig_l(nq2,nom,nta),sig_t(nq2,nom,nta),
     +        sig_lt(nq2,nom,nta),
     +        sig_tt(nq2,nom,nta)
      
      real*8 p_cms,e_cms,p_ara,p_erp,p_mom,p_ene
      real*8 s_cms,w_cms,g_ama,g_abe
      real*8 e_sca,p_sca,c_tet
      real*8 e_bea
      real*8 g_mom
      
      real*8 xm_e,xm_p,xm_n,xm_d,xmd,spn,dpn,xme2,xmn2,xmp2,xmd2,spn2
      real*8 eps,epsp,flux,jcob,sigl,sigt,siglt,sigtt
      real*8 xb_min,xb_max,om_min,om_max
      real*8 cfin,pi
      

      common/fc__some_mas / dpn,xmd,spn2,xmd2,xmn2
      common/fc__sigmas   / sigl,sigt,siglt,sigtt
      
c      common/fc__grid_par / nq2,nom,nta
      common/fc__laget_cnt/ laget_ps_fail,laget_grid_fail
      common/fc__laget_mod/ laget_intp,laget_pwia,laget_fsi,laget_mec
     
      data  xm_e / 510.99890d-03 /,  xm_p / 938.27200d+00 /, 
     +      xm_n / 939.56533d+00 /,  xm_d / 187.56128d+01 /
     
      data  cfin / 137.03599d+00 /
     
      pi = dacos( -1.d+00 )

*---- particle masses and relevant combinations                             

      xmd = xm_d
*-    dimension 1 combinations
      spn = xm_n + xm_p
      dpn = xm_n - xm_p
*-    dimension 2 combinations                                     
      xme2 = xm_e * xm_e                                            
      xmd2 = xm_d * xm_d                                                       
      xmp2 = xm_p * xm_p
      xmn2 = xm_n * xm_n
      spn2 =  spn * spn                                                      

      
*---- initialisations

      xsec  = 0.d+00
      sigl  = 0.d+00
      sigt  = 0.d+00
      siglt = 0.d+00
      sigtt = 0.d+00 
      
*---- beam energy           
                                             
      e_bea = dsqrt( p_bea * p_bea + xme2 )
     
*---- phase space restriction

      xb_min = e_bea * dsqrt(g_mas) + 
     +         p_bea*dsqrt( g_mas + 4.d+00*xme2 )
      xb_min = dsqrt(g_mas) * xb_min / xm_p / 
     +         ( 4.d+00*p_bea*p_bea - g_mas )
      xb_max = xm_d * g_mas / xm_p / ( g_mas + spn2 - xmd2 )
      om_min = 0.5d+00 * g_mas / xm_p / xb_max
      om_max = 0.5d+00 * g_mas / xm_p / xb_min
      if( (g_ene.lt.om_min).or.(g_ene.gt.om_max) ) then
         laget_ps_fail = laget_ps_fail + 1
         write (6,*) 'parameres  : ', xb_min, xb_max, om_min, om_max
         write (6,*) 'Impossible Kinematics '
         write (6,*) 'Beam Energy : ', e_bea
         write (6,*) 'Q2          : ', g_mas
         write (6,*) 'omega       : ', g_ene
         write (6,*) 'theta-p     :', p_tgs
         write (6,*) 'phi-p       : ', p_phi
         go to 1      
      endif

*---- kinematics of the current event  
                                                                              
*-    scattered electron 
      e_sca = e_bea - g_ene                                 ! energy
      p_sca = dsqrt( e_sca * e_sca - xme2 )                 ! momentum
      c_tet = e_bea * e_sca - xme2 - 0.5d+00 * g_mas 
      c_tet = c_tet / p_bea / p_sca                         ! polar angle 
                                                              ! cosinus
*-    virtual photon
      g_mom = dsqrt( g_mas + g_ene * g_ene )                ! momentum
*-    center of mass frame
      s_cms = xmd2 - g_mas + 2.d+00*xm_d*g_ene              ! invariant squared
                                                              ! mass
      w_cms = dsqrt( s_cms )                                ! invariant mass
      g_ama = ( g_ene + xm_d ) / w_cms                      ! gamma
      g_abe = g_mom / w_cms                                 ! gamma * beta
*-    knocked-out proton       
      p_cms = ( s_cms - spn2 ) * ( s_cms - dpn * dpn ) / s_cms
      p_cms = 0.5d+00 * dsqrt( p_cms )                      ! momentum in the 
                                                              ! cms frame
      e_cms = dsqrt( p_cms * p_cms + xmp2 )                 ! energy in the 
                                                              ! cms frame
      p_ara = g_abe * e_cms + g_ama * p_cms * dcos( p_tgs )
      p_erp = p_cms * dsin( p_tgs )
      p_mom = dsqrt( p_ara*p_ara + p_erp*p_erp )            ! momentum in the 
                                                              ! lab frame
      p_ene = g_ama * e_cms + g_abe * p_cms * dcos( p_tgs ) ! energy in the 
                                                              ! lab frame
      
*---- virtual photon polarization

      eps = g_mom * g_mom * (1.d+00 - c_tet) / (1.d+00 + c_tet) / g_mas
      eps = 1.d+00 / ( 1.d+00 + eps + eps ) 
      epsp = dsqrt( 0.5d+00 * g_mas * eps * (1.d+00 + eps) ) / g_ene

*---- virtual photon flux

      flux = 0.5d+00 * e_sca * g_mom / pi / pi / e_bea / g_mas / cfin
      flux = flux / ( 1.d+00 - eps )

*---- calculation of the center of mass to the lab frame jacobian

      jcob = g_ene + xm_d - ( p_ene * g_mom * p_ara / p_mom / p_mom )
      jcob = w_cms * p_mom / p_cms / dabs( jcob )

*---- interpolation of the individual response functions

      if( laget_intp.eq.1 ) then
*-       linear interpolation in angle     
         sigl  = ff__intrpol( sig_l,g_mas,g_ene,p_tgs)
         sigt  = ff__intrpol( sig_t,g_mas,g_ene,p_tgs)
         siglt = ff__intrpol(sig_lt,g_mas,g_ene,p_tgs)
         sigtt = ff__intrpol(sig_tt,g_mas,g_ene,p_tgs)
      else
*-       logarythmic interpolation in recoil momentum
         sigl  = ff__intrlog( sig_l,g_mas,g_ene,p_tgs)
         sigt  = ff__intrlog( sig_t,g_mas,g_ene,p_tgs)
         siglt = ff__intrlog(sig_lt,g_mas,g_ene,p_tgs)
         sigtt = ff__intrlog(sig_tt,g_mas,g_ene,p_tgs)
      endif
*-    safety check to protect against overflow (nan)      
      if((sigl.eq.0.d+00).and.(sigt.eq.0.d+00).and.(siglt.eq.0.d+00)
     +                   .and.(sigtt.eq.0.d+00)) then
         write (6,*) 'all is zero, interpolation problem !'
         go to 1
      endif

*---- differential cross section (ub.mev-1.sr-2)
 
      xsec = flux * jcob 
     +            * ( 
     +              sigt 
     +            + eps  * ( sigl + sigtt * dcos(p_phi+p_phi) )
     +            - epsp *   siglt * dcos(p_phi) 
     +		    )
      
 1    return
      end

c
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c    

      subroutine fs__grid_load(sig_l,sig_t,sig_lt,sig_tt, error)
      
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  author: e. voutier
c  date: october 2005
c  purpose: load the partial cross section data grid
c  
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  load into the common area the  response  functions  of the d(e,e'p) 
c  reaction determined  in 
c  the framework of jean-marc laget's formalism over a grid  sampled  
c  in  q2, omega  and tta_p (the proton angle in the center of mass frame) 
c  in the range:
c
c            0.05 gev2 <= q2 <= 5.00 gev2   step = 0.05 gev2 
c           0.01 gev <= omega <= 5.00 gev   step = 0.01 gev
c                 0 dg <= tta_p <= 180 dg   step = 2 dg 
c
c  the physics options  are  specified  according  to the convention 
c  encoded in the parameters ipwia, ifsi,and imec that are part of the 
c  filename.
c
c     000  neutron contribution only : pwia 
c     001  neutron contribution only : pwia + mec
c     010  neutron contribution only : pwia + fsi
c     011  neutron contribution only : pwia + fsi + mec 
c     100  proton  contribution only : pwia 
c     101  proton  contribution only : pwia + mec
c     110  proton  contribution only : pwia + fsi
c     111  proton  contribution only : pwia + fsi + mec
c     200  neutron + proton contrib. : pwia
c     201  neutron + proton contrib. : pwia + mec 
c     210  neutron + proton contrib. : pwia + fsi 
c     211  neutron + proton contrib. : pwia + fsi + mec
c
c  the physics option parameters are passed via common :
c
c     laget_pwia
c     laget_fsi
c     laget_mec
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - -  -  -  -  -  -  -  -  -
c
c  output variables :
c 
c        sig_l   array of the longitudinal cross sections
c        sig_t   array of the transverse cross sections
c        sig_lt  array of the longitudinal-transverse cross sections
c        sig_tt  array of the transverse-transverse cross sections
c
c  these photoproduction like cross sections are connected to the usual
c  response functions via simple kinematic factors.    
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - -  -  -  -  -  -  -  -  -

      implicit none

      logical save_grid, use_binary
      common/fc__control_c/ save_grid, use_binary

      common/fc__datdir_c/ dat_dir
      common/fc__datdir_i/ n_dat_dir
      common/fc__laget_mod/ laget_intp,laget_pwia,laget_fsi,laget_mec
      
      character*4   q2val,omval
      character*1   pw,fs,mc
      character*300 fname,tmpnam
      character*100 dat_dir
      character*50  subdir

      integer error, i1, i2

      integer*4 n_q2,n_om,n_ta, n_tot
      integer*4 i_q2,i_om, i_count
      integer*4 nchar,n_dat_dir
      integer*4 laget_intp,laget_pwia,laget_fsi,laget_mec

c      parameter ( n_q2 = 100 )
c      parameter ( n_om = 500 )
c      parameter ( n_ta =  91 )

      include 'grid_par.inc'
           
      real*8  sig_l(nq2,nom,nta),  sig_t(nq2,nom,nta),
     +       sig_lt(nq2,nom,nta), sig_tt(nq2,nom,nta)
      real*8 q2_step,om_step,ta_step
      
      real*8 g_mas,g_ene

      real frac
          
      common/fc__grid_deu / q2_step,om_step,ta_step
      
*-    grid parameters

      q2_step = 5.d+06 / dfloat(nq2)
      om_step = 5.d+03 / dfloat(nom)
      ta_step = dacos(-1.d+00) / dfloat(nta-1)
      
      write(pw,'(i1)') laget_pwia
      write(fs,'(i1)') laget_fsi
      write(mc,'(i1)') laget_mec 
            
*---- user selection of the physics grid 
  
*-    no fsi nor mec
      if(laget_fsi.eq.0) then
         subdir = 'deut_laget/pwia/'
*-    with fsi and mec	 
      elseif(laget_mec.eq.1) then
         subdir = 'deut_laget/pful/'
*-    with fsi only	 
      else
         subdir = 'deut_laget/pfsi/'
      endif

      write(6,*) 'interpolate data from : '//subdir

*---- load the physics selected grid	  

      if (use_binary) then
         call fs__nospac(subdir,i1, i2)
         tmpnam = subdir(i1:i2)//'grid.bin'
         call fs__nospac(tmpnam, i1,i2)
         print *, 'using binary grid : '//tmpnam(i1:i2)
         open(10, file = tmpnam(i1:i2), form = 'unformatted', 
     >        status = 'old', err = 998)
         read(10, err = 999) sig_l,sig_t,sig_lt,sig_tt
         error = 0
         close(10)
         return
      endif
 
      n_tot = nq2 * nom
      i_count = 0
      do i_q2=1,nq2
         g_mas = q2_step * dfloat(i_q2)
         write(q2val,'(f4.2)') 1.d-06 * g_mas 
         do i_om=1,nom
            g_ene = om_step * dfloat(i_om)
            write(omval,'(i4)') idint(g_ene)
            if(idint(g_ene).lt.10) then
               omval = '000'//omval(4:4)
            elseif(idint(g_ene).lt.100) then 
               omval = '00'//omval(3:4)
            elseif(idint(g_ene).lt.1000) then 
               omval = '0'//omval(2:4)
            endif
            tmpnam = dat_dir(1:n_dat_dir)//'/'//subdir//
     +               'q2_'//q2val//'_om_'//omval//
     +               '_'//pw//fs//mc//'.dat'
            call fs__squeeze(tmpnam,fname,nchar)
            call fs__read_data
     +           (i_q2,i_om,sig_l,sig_t,sig_lt,sig_tt,fname, error) 
            if (error .ne. 0) then
               return ! error loading data
            endif
            i_count = i_count + 1
            frac = float(i_count)/n_tot*100
            if (mod(frac, 5.) .lt. 1.e-5) then
               print *, frac, ' % read'
            endif
         enddo
      enddo

      print *, 'save_grid = ', save_grid, ' use binary = ', use_binary

      if (save_grid) then
         call fs__nospac(subdir,i1, i2)
         tmpnam = subdir(i1:i2)//'grid.bin'
         call fs__nospac(tmpnam,i1,i2)
         print *, 'saving binary data grid in : '//tmpnam(i1:i2)
         open(10, file = tmpnam(i1:i2), form = 'unformatted')
         write(10) sig_l,sig_t,sig_lt,sig_tt
         close(10)
      endif                                
      return

 998  print *, 'cannot open : '//tmpnam(i1:i2)
      error = -1
      return

 999  print *, 'read error from: '//tmpnam(i1:i2)
      error = -2
      return
      end
      
c
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
c    

      subroutine fs__read_data(ix,iy,xl,xt,xlt,xtt,fname,error)
      
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  author: e. voutier
c  date: october 2005
c  purpose: open and read a basic data file with specified name 
c  
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c  
c  utility  tool for opening the cross section data  base  of jean-marc laget 
c  and load partial cross sections into specified arrays. 
c  
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  input variables : 
c
c        ix      quadrimomentum transfer index         
c        iy      energy transfer index
c        fname   data file name
c
c  output variables :
c 
c        xl      array of the longitudinal cross sections
c        xt      array of the transverse cross sections
c        xlt     array of the longitudinal-transverse cross sections
c        xtt     array of the transverse-transverse cross sections
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

      implicit none

      include 'grid_par.inc'
      
      character*(*) fname
      integer*4 ix,iy,iz
      integer error
        
      real*8 xlt(nq2,nom,nta),xtt(nq2,nom,nta)
      real*8 xl(nq2,nom,nta),xt(nq2,nom,nta)
      
      open(10,file=fname,status='old',err = 998)
      error = 0
      do 100 iz = 1,nta
  100 read(10,'(4(2x,e16.9))',err = 999) 
     +      xl(ix,iy,iz),xt(ix,iy,iz),xlt(ix,iy,iz),xtt(ix,iy,iz) 
c  100 read(10,*,err = 999) 
c     +      xl(ix,iy,iz),xt(ix,iy,iz),xlt(ix,iy,iz),xtt(ix,iy,iz) 
c      close(10) 
  
      return
 998  write (6,*) ' cannot open : ', fname
      error = -1
      return
 999  write (6,*) ' error reading : ', fname
      error = -2
      return
      end

c
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
c    

      double precision function ff__intrpol(xdat,x,y,z)
      
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  author: e. voutier
c  date: october 2005
c  purpose: 3-dimensional linear interpolation 
c  
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c  
c  determine the selected partial cross section (response function) via a 
c  3-dimensional linear interpolation. 
c  
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  input variables : 
c
c        xdat    3-dimensional array of cross section data         
c        x       quadrimomentum transfer
c        y       energy transfer
c        z       proton angle in the enter of mass frame
c
c  output variable :
c
c        ff__intrpol interpolated value of the cross section at (x,y,z)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

      implicit none
      
      include 'grid_par.inc'
        
      integer*4 ix1,iy1,iz1,ix2,iy2,iz2

      real*8 q2_step,om_step,ta_step
      
      real*8 xdat(nq2,nom,nta),x,y,z
      
      real*8 s_111,s_121,s_112,s_122,s_211,s_221,s_212,s_222
      real*8 ax,ay,az,axy,axz,ayz,axyz
      
      common/fc__grid_deu / q2_step,om_step,ta_step
c      common/fc__grid_par / nq2,nom,nta
      
      ff__intrpol = 0.d+00
      
*---- lower array index

      ix1 = idint( x / q2_step )
      iy1 = idint( y / om_step )
      iz1 = idint( z / ta_step ) + 1
      
*---- upper array index

      ix2 = ix1 + 1
      iy2 = iy1 + 1
      iz2 = iz1 + 1
*-    correction for grid boundaries
      if(ix1.eq.nq2) ix2 = ix1
      if(iy1.eq.nom) iy2 = iy1
      if(iz1.eq.nta) iz2 = iz1
      
      if( (ix1.lt.1).or.(ix2.gt.nq2).or.(iy1.lt.1).or.(iy2.gt.nom) )then
         write(6,*) ' cross section grid out of range '
         ff__intrpol = 0.d0
         return
      endif
      
      ax   = (  x - q2_step * dfloat(ix1)   ) / q2_step 
      ay   = (  y - om_step * dfloat(iy1)   ) / om_step 
      az   = (  z - ta_step * dfloat(iz1-1) ) / ta_step
      axy  =   ax * ay
      axz  =   ax * az
      ayz  =   ay * az
      axyz =   ax * ay * az
      
      s_111 = xdat(ix1,iy1,iz1)
      s_121 = xdat(ix1,iy2,iz1)
      s_112 = xdat(ix1,iy1,iz2)
      s_122 = xdat(ix1,iy2,iz2)
      s_211 = xdat(ix2,iy1,iz1)
      s_221 = xdat(ix2,iy2,iz1)
      s_212 = xdat(ix2,iy1,iz2)
      s_222 = xdat(ix2,iy2,iz2)
      
      ff__intrpol =   s_111 * 
     +           ( 1.d+00 - ax - ay - az + axy + axz + ayz - axyz )
     +	      +   s_121 *
     +	        ( ay - axy - ayz + axyz )
     +	      +   s_112 *
     +	        ( az - axz - ayz + axyz )
     +	      +   s_122 *
     +	        ( ayz - axyz )
     +	      +	  s_211 *
     + 	        ( ax - axy - axz + axyz )
     +	      +	  s_221 *
     +	        ( axy - axyz )
     +	      +	  s_212 *
     + 	        ( axz - axyz )	     	
     +	      +	  s_222 *
     +	        ( axyz )
     
      return
      end

c
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
c
 
      double precision function ff__intrlog(xdat,x,y,z)
       
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  author: e. voutier
c  date: october 2005
c  purpose: 3-dimensional interpolation (2 linear + 1 logarythmic)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  determine the selected partial cross section (response function) via a 
c  3-dimensional interpolation assuming a  linear  interpolation in x and y, 
c  and a logarythmic  interpolation  in  f(z) corresponding to the recoil
c  momentum. 
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  input variables : 
c
c        xdat    3-dimensional array of cross section data         
c        x       quadrimomentum transfer
c        y       energy transfer
c        z       proton angle in the center of mass frame
c
c  output variable :
c
c        ff__intrlog interpolated value of the cross section at (x,y,z)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
 
      implicit none
       
      include 'grid_par.inc'
         
      integer*4 ix1,iy1,iz1,ix2,iy2,iz2
 
      real*8 q2_step,om_step,ta_step
       
      real*8 xdat(nq2,nom,nta),x,y,z
      
      real*8 ff__inlgcor
       
      real*8 s_111,s_121,s_112,s_122,s_211,s_221,s_212,s_222
      real*8 q2_1,q2_2,om_1,om_2,pt_1,pt_2
      real*8 s_11,s_12,s_21,s_22
      real*8 ax,ay,axy
        
      common/fc__grid_deu / q2_step,om_step,ta_step
c      common/fc__grid_par / nq2,nom,nta
       
      ff__intrlog = 0.d+00
       
*---- lower array index
 
      ix1 = idint( x / q2_step )
      iy1 = idint( y / om_step )
      iz1 = idint( z / ta_step ) + 1
       
*---- upper array index
 
      ix2 = ix1 + 1
      iy2 = iy1 + 1
      iz2 = iz1 + 1
*-    correction for grid boundaries
      if(ix1.eq.nq2) ix2 = ix1
      if(iy1.eq.nom) iy2 = iy1
      if(iz1.eq.nta) iz2 = iz1
       
      if( (ix1.lt.1).or.(ix2.gt.nq2).or.(iy1.lt.1).or.(iy2.gt.nom) )then
         write(6,*) ' cross section grid out of range '
         ff__intrlog = 0.d0
         return
      endif
       
      q2_1 = q2_step * dfloat(ix1)
      q2_2 = q2_step + dfloat(ix2)
      om_1 = om_step * dfloat(iy1)
      om_2 = om_step * dfloat(iy2)
      pt_1 = ta_step * dfloat(iz1-1)
      pt_2 = ta_step * dfloat(iz2-1)
      
      s_111 = xdat(ix1,iy1,iz1)
      s_121 = xdat(ix1,iy2,iz1)
      s_112 = xdat(ix1,iy1,iz2)
      s_122 = xdat(ix1,iy2,iz2)
      s_211 = xdat(ix2,iy1,iz1)
      s_221 = xdat(ix2,iy2,iz1)
      s_212 = xdat(ix2,iy1,iz2)
      s_222 = xdat(ix2,iy2,iz2)
           
      s_11 = ff__inlgcor(q2_1,om_1,pt_1,z,pt_2,s_111,s_112)
      s_12 = ff__inlgcor(q2_1,om_2,pt_1,z,pt_2,s_121,s_122)
      s_21 = ff__inlgcor(q2_2,om_1,pt_1,z,pt_2,s_211,s_212)
      s_22 = ff__inlgcor(q2_2,om_2,pt_1,z,pt_2,s_221,s_222)
              
      ax   = (  x - q2_1 ) / q2_step
      ay   = (  y - om_1 ) / om_step
      axy  =   ax * ay
       
      ff__intrlog =   s_11 * ( 1.d+00 - ax - ay + axy )
     +        +   s_12 * ( ay - axy )
     +        +   s_21 * ( ax - axy )
     +        +   s_22 * ( axy )
      
      return
      end
 
c
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
c
 
      double precision function ff__inlgcor(x,y,z1,z,z2,f1,f2)
       
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  author: e. voutier
c  date: october 2005
c  purpose: 1-dimensional logarythmic interpolation
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  perform a logarythmic interpolation between the points z1 and z2 for 
c  non-singular f1 value.  singular cases are reduced  to  a  linear 
c  interpolation  and correspond to  a nul f1 (phase space effects at the 
c  kinematic boundaries) and a sign change between f1 and f2.   
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c  input variables : 
c
c        xdat    3-dimensional array of cross section data         
c        x       quadrimomentum transfer
c        y       energy transfer
c        z1      lower bound of the proton angle
c        z       proton angle in the center of mass frame
c        z2      upper bound of the proton angle
c        f1      cross section value at (x,y,z1)
c        f2      cross section value at (x,y,z2)
c
c  output variable :
c
c        ff__inlgcor interpolated value of the function at z
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
      implicit none
      
      real*8 p_recoil
      
      real*8 x,y,z1,z,z2,f1,f2
      
      real*8 dz,f1r,f1s,f2s
      real*8 r,r1,r2,dr

      integer laget_ps_fail,laget_grid_fail

      common/fc__laget_cnt/ laget_ps_fail,laget_grid_fail
      
      ff__inlgcor = f1

      dz = z2 - z1 

*---- grid boundary
      
      if( dz.eq.0.d+00 ) then
         laget_grid_fail = laget_grid_fail + 1
         go to 1
      endif
      
*---- sign extraction

      f1s = 1.d+00
      if( f1.lt.0.d+00 ) f1s = -1.d+00
      f2s = 1.d+00
      if( f2.lt.0.d+00 ) f2s = -1.d+00

*---- recoil neutron momenta
      
      r  = p_recoil(x,y, z)
      r1 = p_recoil(x,y,z1)
      r2 = p_recoil(x,y,z2)
      dr = ( r - r1 ) / ( r2 - r1 )
      if( dr.eq.0.d+00 ) go to 1

*---- quasi-threshold identification

      f1r = 1.d+06
      if( f1.ne.0.d+00 ) f1r = dabs( f2/f1 )

*---- linear interpolation in angle for quasi-threshold effects
      
      if( f1r.gt.1.d+05 ) then
      
         ff__inlgcor = f1 + ( ( f2 - f1 ) * ( z - z1 ) / dz ) 

      elseif( dabs(f1s+f2s).lt.1.d+00 ) then
      
*---- linear interpolation in p_r for oscillating sign

         ff__inlgcor = f1 + ( f2 - f1 ) * dr 
      
      else
	  
*---- logarythmic interpolation for constant sign
*     peu:  check for both signs negative.

         if (f1s+f2s .lt. -1.5) then    ! both signs negative
            ff__inlgcor = -( (-f1)**(1.d+00-dr) ) * ( (-f2)**dr )
         else
            ff__inlgcor = ( f1**(1.d+00-dr) ) * ( f2**dr )
         endif
      
      endif
     
    1 return
      end
       
c
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
c
 
      double precision function p_recoil(g_mas,g_ene,p_tgs)
       
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  author: e. voutier
c  date: october 2005
c  purpose: calculate the 3-momentum of the recoil neutron
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  determine the 3-momentum of the recoiling neutron in the  laboratory frame, 
c  for the current kinematics specified by the virtual photon and the center 
c  of mass angle of the proton.
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
c
c  input variables : 
c
c        g_mas    quadrimomentum transfer         
c        g_ene    energy transfer
c        p_tgs    proton angle in the center of mass frame
c
c  output variable :
c
c        p_recoil corresponding recoil momentum in the lab frame
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
 
      implicit none
      
      real*8 g_mas,g_ene,p_tgs
      
      real*8 dpn,xmd,spn2,xmd2,xmn2
      
      real*8 s_cms,w_cms,r_cms,e_cms
      real*8 g_ama,g_abe,r_par,r_per
      
      common/fc__some_mas / dpn,xmd,spn2,xmd2,xmn2
      
*---- center of mass frame
c      write (6,*) ' dpn = ', dpn,
c     +     ' xmd = ', xmd,
c     +     ' spn2 = ', spn2,
c     +     ' xmd2 = ', xmd2,
c     +     ' xmn2 = ', xmn2

      s_cms = xmd2 - g_mas + 2.d+00*xmd*g_ene                 ! invariant 
                                                                ! squared mass
      w_cms = dsqrt( s_cms )                                  ! invariant mass
      g_ama = ( g_ene + xmd ) / w_cms                         ! gamma
      g_abe = dsqrt( g_mas + g_ene * g_ene ) / w_cms          ! gamma * beta

*---- recoil neutron in the center of mass frame

      r_cms = ( s_cms - spn2 ) * ( s_cms - dpn * dpn ) / s_cms
      r_cms = 0.5d+00 * dsqrt( r_cms )                        ! momentum in the
                                                                ! cms frame
      e_cms = dsqrt( r_cms * r_cms + xmn2 )                   ! energy in the 
                                                                ! cms frame
      
*---- recoil neutron in the lab frame
      
      r_par =   g_abe * e_cms - g_ama * r_cms * dcos( p_tgs ) ! parallel mom.
                                                                ! component
      r_per = - r_cms * dsin( p_tgs )                         ! perp. mom.
                                                                ! component 
      
c      write (6,*) 'p_recoil :', ' gamma_cm = ', g_ama, 
c     +            ' gamma*beta = ', g_abe, 
c     +            ' beta = ', g_abe/g_ama, 
c     +            ' e_cms = ', e_cms,
c     +            ' r_per = ', r_per,
c     +            ' r_par = ', r_par,
c     +            ' r_par_cm = ', r_cms*dcos(p_tgs)
      p_recoil = dsqrt( r_par*r_par + r_per*r_per )           ! recoil momentum
                                                                ! in lab frame
      
      return
      end
      
c
c-----------------------------------------------------------------------------
c  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  
c-----------------------------------------------------------------------------
c
c------------------------------------------------------------------------------
c       extracts non-blank characters of any character string.
c------------------------------------------------------------------------------
c
      subroutine fs__squeeze(charin,charout,nchar)
      implicit none
c
      integer nchar,n,j
      character*300 charin,charout
      n = 0
      do j=1,300
        if(charin(j:j).ne.' ') then
          n = n+1
          charout(n:n) = charin(j:j)
        endif
      enddo
c
      do j=n+1,300
         charout(j:j) = ' '
      enddo
c
      nchar = n
      return
      end
c
      subroutine fs__nospac(ch,i1,i2)
c     find first and last non-space characters in an arbitrary length
c     string.  if string is only spaces, i2 will be 1 less than i1
      character ch*(*)
      i2=len(ch)
      i1=1
6000  if (i1.lt.i2.and.ch(i1:i1).eq.' ') then
          i1=i1+1
      goto 6000
      endif
6005  if (i2.ge.i1.and.ch(i2:i2).eq.' ') then
          i2=i2-1
      goto 6005
      endif
      return
      end
