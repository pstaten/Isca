module damping_driver_mod

    !-----------------------------------------------------------------------
    !
    !       This module controls four functions:
    !
    !   (1) rayleigh friction applied to momentum fields at levels
    !       1 to kbot (i.e., momentum is damped toward zero).
    !
    !   (2) mountain gravity wave drag module may be called
    !
    !   (3) Alexander-Dunkerton gravity wave drag may be called
    !
    !   (4) Garner topo_drag module may be called
    !mj
    !   (5) Time independent "gravity wave" drag may be called
    !
    !-----------------------------------------------------------------------
    
     use      mg_drag_mod, only:  mg_drag, mg_drag_init, mg_drag_end
     use      cg_drag_mod, only:  cg_drag_init, cg_drag_calc, cg_drag_end
     use    topo_drag_mod, only:  topo_drag_init, topo_drag, topo_drag_end
     use          fms_mod, only:  file_exist, mpp_pe, mpp_root_pe, stdlog, &
                                  write_version_number, &
                                  open_namelist_file, error_mesg, &
                                  check_nml_error,                   &
                                  FATAL, close_file
     use diag_manager_mod, only:  register_diag_field,  &
                                  register_static_field, send_data
     use time_manager_mod, only:  time_type,get_time,length_of_year !mj
     use    constants_mod, only:  cp_air, grav, PI
    
     implicit none
     private
    
     public   damping_driver, damping_driver_init, damping_driver_end
    
    !-----------------------------------------------------------------------
    !---------------------- namelist ---------------------------------------
    
       real     :: trayfric = 0.
    ! mj pk02-like sponge   integer  :: nlev_rayfric = 1
       integer  :: nlev_rayfric
       real :: sponge_pbottom = 50. ! [Pa]
      logical  :: do_mg_drag = .false.
    !epg: Use cg_drag.f90, GFDL's version of the Alexander and Dunkerton 1999 
    !     Non-orographic gravity wave parameterization, updated as for Cohen et al. 2013
    ! mj actively choose rayleigh friction
       logical  :: do_rayleigh = .false.
       logical  :: do_cg_drag = .false.
       logical  :: do_topo_drag = .false.
       logical  :: do_const_drag = .false.
       real     :: const_drag_amp = 3.e-04
       real     :: const_drag_off = 0.
       logical  :: do_conserve_energy = .false.
    ! also add the fake_qbo subroutine
       logical  :: do_sin_qbo = .true.
       logical  :: do_ewa_htg = .false.

    ! namelist parameters for ewa_heating
       real     :: h_amp = 0.1/86400.     ! Kelvin per day
       real     :: p_s = 100000.         ! Pascals
       real     :: p_center = 5000.      ! Pascals
       real     :: p_width = 0.0175      ! Pascals
       real     :: lat_width = 0.4       ! Radians

    ! namelist parameters for fake_qbo
       real     :: qbo_amp = 20.               ! m/s

    
       namelist /damping_driver_nml/  trayfric,  &
                                      do_rayleigh, sponge_pbottom,  & ! mj
                                      do_cg_drag, do_topo_drag, &
                                      do_mg_drag, do_conserve_energy, &
                                      do_const_drag, const_drag_amp,const_drag_off, &
                                      do_sin_qbo, do_ewa_htg, &
                                      h_amp, p_s, p_center, p_width, lat_width, A
                                          !mj
    
    !
    !   trayfric = damping time in seconds for rayleigh damping momentum
    !              in the top nlev_rayfric layers (if trayfric < 0 then time
    !              in days)
    !                 [real, default: trayfric=0.]
    !
    !   nlev_rayfric = number of levels at the top of the model where
    !                  rayleigh friction of momentum is performed, if
    !                  trayfric=0. then nlev_rayfric has no effect
    !                    [integer, default: nlev_rayfric=1]
    !
    !-----------------------------------------------------------------------
    !----- id numbers for diagnostic fields -----
    
    integer :: id_udt_rdamp,  id_vdt_rdamp,   &
               id_udt_gwd,    id_vdt_gwd,     &
                              id_sgsmtn,      &
               id_udt_cgwd,   id_taus,        &
               id_udt_cnstd                 !mj
    
    integer :: id_tdt_diss_rdamp,  id_diss_heat_rdamp, &
               id_tdt_diss_gwd,    id_diss_heat_gwd
    
    integer :: id_udt_topo,   id_vdt_topo,   id_taubx,  id_tauby
    
    !----- missing value for all fields ------
    
    real :: missing_value = -999.
    
    character(len=7) :: mod_name = 'damping'
    
    !-----------------------------------------------------------------------
    !mj actively choose rayleigh - is now in namelist
    ! logical :: do_rayleigh
    
     real, parameter ::  daypsec=1./86400.
     logical :: module_is_initialized =.false.
    
     real :: rfactr
    
    !   note:  
    !     rfactr = coeff. for damping momentum at the top level
    
     character(len=128) :: version = '$Id: damping_driver.f90,v 10.0 2003/10/24 22:00:25 fms Exp $'
     character(len=128) :: tagname = '$Name: lima $'
    
    !mj cg_drag alarm
     integer :: Time_lastcall,dt_integer,days,seconds
    !-----------------------------------------------------------------------
    
    contains
    
    !#######################################################################
    
     subroutine damping_driver (is, js, lat, Time, delt, pfull, phalf, zfull, zhalf, &
                                u, v, t, q, r,  udt, vdt, tdt, qdt, rdt,  &
    !                                   mask, kbot)
                                z_pbl,  mask, kbot)
     
    !-----------------------------------------------------------------------
     integer,         intent(in)                :: is, js
     real, dimension(:,:), intent(in)           :: lat
     type(time_type), intent(in)                :: Time
     real,            intent(in)                :: delt
     real,    intent(in),    dimension(:,:,:)   :: pfull, phalf, &
                                                   zfull, zhalf, &
                                                   u, v, t, q
     real,    intent(in),    dimension(:,:,:,:) :: r
     real,    intent(inout), dimension(:,:,:)   :: udt,vdt,tdt,qdt
     real,    intent(inout), dimension(:,:,:,:) :: rdt
     real, dimension(:,:), intent(in)           :: z_pbl
     real,    intent(in),    dimension(:,:,:), optional :: mask
     integer, intent(in),    dimension(:,:),   optional :: kbot
     integer :: nlon, nlat, nlev 

    
    !-----------------------------------------------------------------------
     real, dimension(size(udt,1),size(udt,2))             :: diag2
     real, dimension(size(udt,1),size(udt,2))             :: taubx, tauby
     real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: taus
     real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: utnd, vtnd, &
                                                             ttnd, pmass, &
                                                             p2
     logical :: used
    
     real, dimension(size(udt,1),size(udt,2),size(udt,3)+1) :: p_pass, &
                                                                t_pass
     integer :: k, j, i, locmax(3)
     real :: a,b
    !-----------------------------------------------------------------------
    !mj constant drag TOA
     real :: minp,cosday
     integer :: seconds,days,daysperyear
    !-----------------------------------------------------------------------
    
       if (.not.module_is_initialized) call error_mesg ('damping_driver',  &
                         'damping_driver_init must be called first', FATAL)
    
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !----------------- r a y l e i g h   d a m p i n g ---------------------
    !-----------------------------------------------------------------------
       if (do_rayleigh) then
    
    ! mj pk02-like sponge       p2 = pfull * pfull
    !       call rayleigh (delt, p2, u, v, utnd, vtnd, ttnd)
           call rayleigh (delt, pfull, u, v, utnd, vtnd, ttnd)
           udt = udt + utnd
           vdt = vdt + vtnd
           tdt = tdt + ttnd
    
    !----- diagnostics -----
    
           if ( id_udt_rdamp > 0 ) then
    !            used = send_data ( id_udt_rdamp, utnd, Time, is, js, 1, &
    !                               rmask=mask )
                used = send_data ( id_udt_rdamp, utnd, Time, rmask=mask )
           endif
    
           if ( id_vdt_rdamp > 0 ) then
    !            used = send_data ( id_vdt_rdamp, vtnd, Time, is, js, 1, &
    !                               rmask=mask )
                used = send_data ( id_vdt_rdamp, vtnd, Time, rmask=mask )
           endif
    
           if ( id_tdt_diss_rdamp > 0 ) then
    !            used = send_data ( id_tdt_diss_rdamp, ttnd, Time, is, js, 1, &
    !                               rmask=mask )
                used = send_data ( id_tdt_diss_rdamp, ttnd, Time, rmask=mask )
           endif
    
           if ( id_diss_heat_rdamp > 0 ) then
                do k = 1,size(u,3)
                  pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
                enddo
                diag2 = cp_air/grav * sum(ttnd*pmass,3)
    !            used = send_data ( id_diss_heat_rdamp, diag2, Time, is, js )
                used = send_data ( id_diss_heat_rdamp, diag2, Time)
           endif
    
       endif
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !--------- m t n   g r a v i t y   w a v e   d r a g -------------------
    !-----------------------------------------------------------------------
       if (do_mg_drag) then
    
           call mg_drag (is, js, delt, u, v, t, pfull, phalf, zfull, zhalf,  &
                         utnd, vtnd, ttnd, taubx,tauby,taus,        kbot)
           udt = udt + utnd
           vdt = vdt + vtnd
           tdt = tdt + ttnd
    
    !----- diagnostics -----
    
           if ( id_udt_gwd > 0 ) then
    !            used = send_data ( id_udt_gwd, utnd, Time, is, js, 1, &
    !                               rmask=mask )
                used = send_data ( id_udt_gwd, utnd, Time, rmask=mask )
           endif
    
           if ( id_vdt_gwd > 0 ) then
    !            used = send_data ( id_vdt_gwd, vtnd, Time, is, js, 1, &
    !                               rmask=mask )
                used = send_data ( id_vdt_gwd, vtnd, Time, rmask=mask )
           endif
    
           if ( id_taubx > 0 ) then
    !            used = send_data ( id_taubx, taubx, Time, is, js )
                used = send_data ( id_taubx, taubx, Time)
           endif
    
           if ( id_tauby > 0 ) then
    !            used = send_data ( id_tauby, tauby, Time, is, js )
                used = send_data ( id_tauby, tauby, Time)
           endif
    
           if ( id_taus > 0 ) then
    !           used = send_data ( id_taus, taus, Time, is, js, 1, &
    !                              rmask=mask )
               used = send_data ( id_taus, taus, Time, rmask=mask )
           endif
    
           if ( id_tdt_diss_gwd > 0 ) then
    !            used = send_data ( id_tdt_diss_gwd, ttnd, Time, is, js, 1, &
    !                               rmask=mask )
                used = send_data ( id_tdt_diss_gwd, ttnd, Time, rmask=mask )
           endif
    
           if ( id_diss_heat_gwd > 0 ) then
                do k = 1,size(u,3)
                  pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
                enddo
                diag2 = cp_air/grav * sum(ttnd*pmass,3)
    !            used = send_data ( id_diss_heat_gwd, diag2, Time, is, js )
                used = send_data ( id_diss_heat_gwd, diag2, Time )
           endif
    
       endif
    
    !   Alexander-Dunkerton gravity wave drag
    
       if (do_cg_drag) then
    !mj updating call to riga version of cg_drag
          !call cg_drag_calc (is, js, lat, pfull, zfull, t, u, Time,    &
          !                  delt, utnd)
          call cg_drag_calc (is, js, lat, pfull, zfull, t, u, v, Time, delt, utnd, vtnd)
         udt =  udt + utnd
         vdt =  vdt + vtnd !mj
    
    !----- diagnostics -----
    
         if ( id_udt_cgwd > 0 ) then
    !        used = send_data ( id_udt_cgwd, utnd, Time, is, js, 1, &
    !                          rmask=mask )
            used = send_data ( id_udt_cgwd, utnd, Time, rmask=mask )
         endif
     
       endif
    
    ! constant drag, modeled on Alexander-Dunkerton winter average
       if (do_const_drag) then
          ! get time of the year for seasonal cycle
          call get_time(length_of_year(),seconds,daysperyear)
          call get_time(Time,seconds,days)
          cosday = cos(2*PI*days/daysperyear)
          utnd = 0.
          minp = log(minval(pfull*0.01)) - 1.
          where( pfull*0.01 < exp(1.) )
             ! vertical: linear in ln(p)
             utnd = -const_drag_amp*((log(pfull*0.01) - 1.)/minp)**1.
          endwhere
          ! latitudinal: 3rd order polynomial, and cosine in time
          do k=1,size(utnd,3)
             where( pfull(:,:,k)*0.01 < exp(1.) )
                utnd(:,:,k) = utnd(:,:,k)*sign(1.,lat)*cosday &
                     *( -1.65*abs(lat)**3 +2.5*lat**2 +0.17*abs(lat) +const_drag_off )
             end where
          enddo
          udt = udt + utnd 
    
    !----- diagnostics -----
    
         if ( id_udt_cnstd > 0 ) then
    !        used = send_data ( id_udt_cnstd, utnd, Time, is, js, 1, &
    !                          rmask=mask )
            used = send_data ( id_udt_cnstd, utnd, Time, rmask=mask )
         endif
       endif
   
    ! -------call for the use of do_sin_qbo-----------------
       if (do_sin_qbo) then
         print *, 'do_sin_qbo is set'
         nlon = size(pfull,1)
         nlat = size(pfull,2)
         nlev = size(pfull,3)
         call get_time(length_of_year(),seconds,daysperyear) 
         call get_time(Time,seconds,days)
         call fake_qbo(pfull, zfull, lat, u, utnd, seconds, days, daysperyear, nlon, nlat, nlev)    
         udt = udt + utnd
       else
         print *, 'do_sin_qbo is not set'
       end if

       ! -------call for the use of do_ewa_htg -----------------
    if (do_ewa_htg) then
            print *, 'do_ewa_htg is set'
         call ewa_heating(lat,pfull, ttnd)    
         tdt = tdt + ttnd
       else
               print *, 'do_ewa_htg is not set'
       end if
    ! ---------------------------


    !-----------------------------------------------------------------------
    !---------topographic   w a v e   d r a g -------------------
    !-----------------------------------------------------------------------
       if (do_topo_drag) then
    
        call topo_drag ( is, js, u, v, t, pfull, phalf, zfull, zhalf,  &
    !               taubx, tauby, utnd, vtnd,taus)
                    z_pbl, taubx, tauby, utnd, vtnd,taus)
    
         b = maxval(abs(utnd))
         locmax = maxloc(abs(utnd))
    
    
         udt = udt + utnd
         vdt = vdt + vtnd
    
    
    !----- diagnostics -----
    
        if ( id_udt_topo > 0 ) then
    !       used = send_data ( id_udt_topo, utnd, Time, is, js, 1, &
    !                          rmask=mask )
           used = send_data ( id_udt_topo, utnd, Time, rmask=mask )
        endif
    
        if ( id_vdt_topo > 0 ) then
    !         used = send_data ( id_vdt_topo, vtnd, Time, is, js, 1, &
    !                         rmask=mask )
             used = send_data ( id_vdt_topo, vtnd, Time, rmask=mask )
       endif
    
         if ( id_taubx > 0 ) then
    !       used = send_data ( id_taubx, taubx, Time, is, js )
           used = send_data ( id_taubx, taubx, Time)
         endif
    
         if ( id_tauby > 0 ) then
    !        used = send_data ( id_tauby, tauby, Time, is, js )
            used = send_data ( id_tauby, tauby, Time)
          endif
    
         if ( id_taus > 0 ) then
    !      used = send_data ( id_taus, taus, Time, is, js, 1, &
    !                          rmask=mask )
          used = send_data ( id_taus, taus, Time, rmask=mask )
         endif
    
    
    
     endif
    
    !-----------------------------------------------------------------------
    
     end subroutine damping_driver
    
    !#######################################################################
    
     subroutine damping_driver_init ( lonb, latb, pref, axes, Time, sgsmtn)
    
     real,            intent(in) :: lonb(:), latb(:), pref(:)
     integer,         intent(in) :: axes(4)
     type(time_type), intent(in) :: Time
     real, dimension(:,:), intent(out) :: sgsmtn
    !-----------------------------------------------------------------------
    !     lonb  = longitude in radians of the grid box edges
    !     latb  = latitude  in radians of the grid box edges
    !     axes  = axis indices, (/x,y,pf,ph/)
    !               (returned from diag axis manager)
    !     Time  = current time (time_type)
    !     sgsmtn = subgrid scale topography variance
    !-----------------------------------------------------------------------
     integer :: unit, ierr, io
     logical :: used
    !mj
     integer :: raylev(1)
    !-----------------------------------------------------------------------
    !----------------- namelist (read & write) -----------------------------
    
       if (file_exist('input.nml')) then
          unit = open_namelist_file ()
          ierr=1; do while (ierr /= 0)
             read  (unit, nml=damping_driver_nml, iostat=io, end=10)
             ierr = check_nml_error (io, 'damping_driver_nml')
          enddo
     10   call close_file (unit)
       endif
    
       call write_version_number(version, tagname)
       if(mpp_pe() == mpp_root_pe() ) then
            write (stdlog(),nml=damping_driver_nml)
       endif
    
    !-----------------------------------------------------------------------
    !--------- rayleigh friction ----------
    
    !mj actively choose rayleigh friction
    !   do_rayleigh=.false.
    
    !   if (abs(trayfric) > 0.0001 .and. nlev_rayfric > 0) then
       if ( do_rayleigh ) then
    ! mj automatically determine nlev_rayfric
          raylev = minloc(abs(pref(:)-2*sponge_pbottom))
          nlev_rayfric = raylev(1)
          if (trayfric > 0.0) then
             rfactr=(1./trayfric)
          else
             rfactr=(1./abs(trayfric))*daypsec
          endif
       endif
    !      do_rayleigh=.true.
    !   else
    !      rfactr=0.0
    !   endif
    
    !-----------------------------------------------------------------------
    !----- mountain gravity wave drag -----
    
       if (do_mg_drag) call mg_drag_init (lonb, latb, sgsmtn)
    
    !--------------------------------------------------------------------
    !----- Alexander-Dunkerton gravity wave drag -----
     
       if (do_cg_drag)  then
         call cg_drag_init (lonb, latb, pref, Time=Time, axes=axes)
       endif
    
    !-----------------------------------------------------------------------
    !----- initialize diagnostic fields -----
    
    if (do_rayleigh) then
    
       id_udt_rdamp = &
       register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time,       &
                           'u wind tendency for Rayleigh damping', 'm/s2', &
                            missing_value=missing_value               )
    
       id_vdt_rdamp = &
       register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time,       &
                           'v wind tendency for Rayleigh damping', 'm/s2', &
                            missing_value=missing_value               )
    
       id_tdt_diss_rdamp = &
       register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), Time,  &
                          'Dissipative heating from Rayleigh damping',&
                                 'deg_k/s', missing_value=missing_value   )
           
       id_diss_heat_rdamp = &
       register_diag_field ( mod_name, 'diss_heat_rdamp', axes(1:2), Time,   &
                    'Integrated dissipative heating from Rayleigh damping',&
                      'W/m2' )
    endif
    
    if (do_mg_drag) then
    
     ! register and send static field
       id_sgsmtn = &
       register_static_field ( mod_name, 'sgsmtn', axes(1:2), &
                   'sub-grid scale topography for gravity wave drag', 'm')
       if (id_sgsmtn > 0) used = send_data (id_sgsmtn, sgsmtn, Time)
    
     ! register non-static field
       id_udt_gwd = &
       register_diag_field ( mod_name, 'udt_gwd', axes(1:3), Time,        &
                         'u wind tendency for gravity wave drag', 'm/s2', &
                            missing_value=missing_value               )
    
       id_vdt_gwd = &
       register_diag_field ( mod_name, 'vdt_gwd', axes(1:3), Time,        &
                         'v wind tendency for gravity wave drag', 'm/s2', &
                            missing_value=missing_value               )
    
       id_taubx = &
       register_diag_field ( mod_name, 'taubx', axes(1:2), Time,        &
                             'x base flux for grav wave drag', 'kg/m/s2', &
                             missing_value=missing_value               )
    
       id_tauby = &
       register_diag_field ( mod_name, 'tauby', axes(1:2), Time,        &
                             'y base flux for grav wave drag', 'kg/m/s2', &
                             missing_value=missing_value )
    
       id_taus = &
       register_diag_field ( mod_name, 'taus', axes(1:3), Time,        &
                           'saturation flux for gravity wave drag', 'kg/m/s2', &
                          missing_value=missing_value               )
    
       id_tdt_diss_gwd = &
       register_diag_field ( mod_name, 'tdt_diss_gwd', axes(1:3), Time,    &
                              'Dissipative heating from gravity wave drag',&
                                  'deg_k/s', missing_value=missing_value   )
           
       id_diss_heat_gwd = &
       register_diag_field ( mod_name, 'diss_heat_gwd', axes(1:2), Time,      &
                    'Integrated dissipative heating from gravity wave drag',&
                                     'W/m2' )
    endif
    
       if (do_cg_drag) then
    
        id_udt_cgwd = &
        register_diag_field ( mod_name, 'udt_cgwd', axes(1:3), Time,        &
                     'u wind tendency for cg gravity wave drag', 'm/s2', &
                          missing_value=missing_value               )
       endif
    
       if (do_const_drag) then
    
        id_udt_cnstd = &
        register_diag_field ( mod_name, 'udt_cnstd', axes(1:3), Time,        &
                     'u wind tendency for constant drag', 'm/s2', &
                          missing_value=missing_value               )
       endif
          
    
    !-----------------------------------------------------------------------
    !----- topo wave drag -----
    
    
    
      if (do_topo_drag) then
              call topo_drag_init (lonb, latb, ierr)
              sgsmtn(:,:) = -99999.
      endif
    
    
    
      if (do_topo_drag) then
    
       id_udt_topo = &
       register_diag_field ( mod_name, 'udt_topo', axes(1:3), Time,        &
                           'u wind tendency for topo wave drag', 'm/s2', &
                            missing_value=missing_value               )
    
      id_vdt_topo = &
       register_diag_field ( mod_name, 'vdt_topo', axes(1:3), Time,        &
                           'v wind tendency for topo wave drag', 'm/s2', &
                             missing_value=missing_value               )
    
       id_taubx = &
       register_diag_field ( mod_name, 'taubx', axes(1:2), Time,        &
                         'x base flux for topo wave drag', 'kg/m/s2', &
                            missing_value=missing_value               )
    
        id_tauby = &
        register_diag_field ( mod_name, 'tauby', axes(1:2), Time,        &
                        'y base flux for topo wave drag', 'kg/m/s2', &
                          missing_value=missing_value )
    
        id_taus = &
       register_diag_field ( mod_name, 'taus', axes(1:3), Time,        &
                      'saturation flux for topo wave drag', 'kg/m/s2', &
                         missing_value=missing_value               )
    
     endif
    
    
    
    !-----------------------------------------------------------------------
    
       module_is_initialized =.true.
    
    !******************** end of initialization ****************************
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    
     end subroutine damping_driver_init
    
    !#######################################################################
    
     subroutine damping_driver_end
    
         if (do_mg_drag) call mg_drag_end
         if (do_cg_drag)   call cg_drag_end
         if (do_topo_drag) call topo_drag_end
    
         module_is_initialized =.false.
    
    
     end subroutine damping_driver_end
    
    !#######################################################################
    
     subroutine rayleigh (dt, p2, u, v, udt, vdt, tdt)
    
      real,    intent(in)                      :: dt
      real,    intent(in),  dimension(:,:,:)   :: p2, u, v
      real,    intent(out), dimension(:,:,:)   :: udt, vdt, tdt
    
      real, dimension(size(u,1),size(u,2)) :: fact
      integer :: k
    !-----------------------------------------------------------------------
    !--------------rayleigh damping of momentum (to zero)-------------------
    
       udt = 0.
       vdt = 0.
       do k = 1, nlev_rayfric
          where ( p2(:,:,k) < sponge_pbottom ) !mj note: p2==pfull now
    !         fact(:,:) = rfactr*(1.+(p2(:,:,1)-p2(:,:,k))/(p2(:,:,1)+p2(:,:,k)))
             fact(:,:) = rfactr*(sponge_pbottom-p2(:,:,k))**2/(sponge_pbottom)**2
             udt(:,:,k) = -u(:,:,k)*fact(:,:)
             vdt(:,:,k) = -v(:,:,k)*fact(:,:)
          endwhere
       enddo
    !   do k = nlev_rayfric+1, size(u,3)
    !     udt(:,:,k) = 0.0
    !     vdt(:,:,k) = 0.0
    !   enddo
    
    !  total energy conservation
    !  compute temperature change loss due to ke dissipation
    
       tdt = 0. !mj
       if (do_conserve_energy) then
           do k = 1, nlev_rayfric
              tdt(:,:,k) = -((u(:,:,k)+.5*dt*udt(:,:,k))*udt(:,:,k) +  &
                             (v(:,:,k)+.5*dt*vdt(:,:,k))*vdt(:,:,k)) / cp_air
           enddo
    !       do k = nlev_rayfric+1, size(u,3)
    !          tdt(:,:,k) = 0.0
    !       enddo
    !   else
    !       tdt = 0.0
       endif
    !-----------------------------------------------------------------------
    
     end subroutine rayleigh
    
    !#######################################################################

     subroutine ewa_heating(lat,p_full,ttnd)
      real, dimension(:,:)  ,intent(in)    :: lat
      real, dimension(:,:,:),intent(in)    :: p_full
      real, dimension(:,:,:),intent(inout) :: ttnd
      ! local variables
      integer :: i,j,k
      real    :: p_factor,lat_factor

      ! This is an implimentation of the heating profile specified in Ewa Bednarz' work.
      ! Q(\textup{lat}, p)=Q_{max} \cdot \exp{\left [-\left ( \frac{\textup{lat}^2}{2\times0.4^2} +\frac{(p/p_s-50/p_s)^2}{2\times0.0175^2}\right )\right ]}
      ! where Q_{max}=1/8640 K/s (that's 10 K per day, which is big)
      ! and p_s=1000 hPa




      ttnd = 0.
      do k = 1,size(p_full,3)
        do j = 1,size(lat,2)
          do i = 1,size(lat,1)
            ! vertical component
            p_factor = exp( -(p_full(i,j,k)/p_s-p_center/p_s)**2/(2*p_width**2))

            ! latitudinal component
            lat_factor = exp( -(lat(i,j))**2/(2*(lat_width)**2) )

            ! everything together
            ttnd(i,j,k) = h_amp*p_factor*lat_factor
          enddo
        enddo
      enddo
  
     end subroutine ewa_heating

    
    !#######################################################################

     subroutine fake_qbo (pfull, zfull, lat, u, utnd, seconds, days, daysperyear, nlon, nlat, nlev)
        implicit none
        integer, intent(in)                      :: seconds, days, daysperyear
        real,    intent(in),  dimension(nlon,nlat)     :: lat
        real,    intent(in),  dimension(nlon,nlat,nlev)   :: pfull, zfull, u
        real,    intent(out), dimension(nlon,nlat,nlev)   :: utnd
        integer, intent(in) :: nlon, nlat, nlev

        ! qbo: QBO profile (m/s)
        ! alt_profile: altitude profile (nlon x nlat x nlev)
        ! relax_rate: relaxation rate (inverse of nudging timescale) (s^-1)
        real, dimension(nlon,nlat,nlev) :: qbo, alt_profile, relax_rate

        ! lat_profile: latitude profile (nlon x nlat)
        real, dimension(nlon,nlat) :: lat_profile

        ! ilev: loop index for vertical levels
        integer :: ilev

        ! m: QBO vertical wavenumber (m^-1)
        ! nu: QBO frequency (s^-1)
        ! p_b: bottom level of relaxation profile (Pa)
        ! del_p_b: transition width in TANH for bottom level (Pa)
        ! p_t: top level of relaxation profile (Pa)
        ! del_p_t: transition width in TANH for top level (Pa)
        ! lat_n: northern latitude limit of the relaxation profile (radians)
        ! lat_s: southern latitude limit of the relaxation profile (radians)
        ! del_lat: transition width in TANH for latitude profile (radians)
        ! alpha: relaxation rate (s^-1)
        real :: m, nu, p_b, del_p_b, p_t, del_p_t, lat_n, lat_s, del_lat, alpha
        !-----------------------------------------------------------------------

        print *, 'subroutine fake_qbo'

      ! QBO function
      ! qbo_amp is now a namelist parameter
        nu = 2. * acos(-1.0) / ((28.0/12.0) * daysperyear/daypsec) ! 28-month period
        m = -2. * acos(-1.0) / 40000 ! 40 km vertical wavelength
        qbo = qbo_amp * COS(m * zfull - nu * (days*24*3600 + seconds))

      ! Relaxation rate (inverse of nudging timescale)
        alpha = (1./(5 * 24 * 3600)) ! 5-day relaxation timescale

      ! Altitude profile (nlon x nlat x nlev)
        p_b = 8000 ! 80 hPa
        del_p_b = 500 ! 5 hPa
        p_t = 800 ! 8 hPa
        del_p_t = 50 ! 0.5 hPa
        alt_profile = (1-TANH((pfull-p_b)/del_p_b)) * (1+TANH((pfull-p_t)/del_p_t))

      ! Latitude profile (nlon x nlat)
        lat_n = (acos(-1.)*12)/180. ! 12°N
        lat_s = (acos(-1.)*(-12.))/180. ! 12°S
        del_lat = (acos(-1.)*2)/180. ! 2°
        lat_profile = (1+TANH((lat-lat_s)/del_lat)) * (1-TANH((lat-lat_n)/del_lat))

      ! element-wise relaxation rate
        do ilev = 1, nlev,1
          relax_rate(:,:,ilev) = alpha * alt_profile(:,:,ilev) * lat_profile / 16
        end do

      ! Calculate the actual u-tendency from the QBO nudging
        utnd = (qbo - u) * (1.0-EXP(-relax_rate))

    end subroutine fake_qbo

    !#######################################################################
    
    end module damping_driver_mod

