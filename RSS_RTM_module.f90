! Remote Sensing System
! Santa Rosa, California
! 
! Contact 
! Thomas Meissner
! meissner@remss.com
!
! September 14, 2022
!
! RSS_RTM_2012
! published on GitHub under MIT license 
!
! If the RTM plays an essential role in the generation of new research results, please acknowledge/credit RSS and let us know.
!
!
!
! FORTRAN 90 module for computing ocean surface emissivity at microwave frequencies
! Citation:
! T. Meissner and F.  Wentz, RSS Ocean Surface Emissivity Model, doi: 10.21982/M84S3Q.
!
!
! References:
! [MW 2004]:   T. Meissner and F. Wentz, 
!              "The Complex Dielectric Constant of Pure and Sea Water from Microwave Satellite Observations", 
!              IEEE Transactions on Geoscience and Remote Sensing, vol 42 (9), 1836-1849, doi: 10.1109/TGRS.2004.831888, 2004.               
!
! [MW 2012]:   T. Meissner, and F. Wentz, 
!              "The Emissivity of the Ocean Surface between 6 - 90 GHz over a Large Range of 
!              Wind Speeds and Earth Incidence Angles", 
!              IEEE Transactions on Geoscience and Remote Sensing, vol 50 (8), 3004-3026, doi: 10.1109/TGRS.2011.2179662, 2012.     
!
! 
! [MWR 2014]:  T. Meissner, F. Wentz, and L. Ricciardulli, 
!              "The Emission and Scattering of L-band Microwave Radiation 
!              from Rough Ocean Surfaces and Wind speed Measurements from Aquarius", 
!              Journal of Geophysical Research: Oceans, 119, doi:10.1002/2014JC009837, 2014.    
!
! [MWD 2012]:  T. Meissner, F. Wentz and D. Draper, 
!              "GMI Calibration Algorithm and Analysis Theoretical Basis Document", 
!              report number 041912, Version-G, Remote Sensing Systems, Santa Rosa, CA, 124 pp.
!              doi: 10.56236/RSS-au,    
!              https://doi.org/10.56236/RSS-au 
!
!
! Release Notes:
!
! The following changes/updates from [MW 2012] are contained in the code.
! 1. dielectric_meissner_wentz: Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
! 2. dielectric_meissner_wentz: Changed SST behavior of coefficient b2 from:
!     b2 = 1.0 + s*(z(10) + z(11)*sst) to
!     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
! 3. Introduced SST dependence for wind direction signal similar than eq. (15) in [MW 2012].  
! 
!
!
!
! Routines:

! 1.  find_surface_tb                Master wrapper routine that caclulates all components of the ocean surface emissivity. Output user selected.

! 2.  dielectric_meissner_wentz      Dielectric model of sea and pure water [MW 2004], with minor updates in [MW 2012].

! 3.  fdem0_meissner_wentz           Caclulates emissivity of specular surface [MW 2004]

! 4.  fd_emiss                       Calculates emissivity of specular surface (v/h), isotropic wind induced emissivity (v/h) and wind direction signal (v/h/S3/S4) 
!                                    [MW 2012], sections IV + VI.
!
! 5.  fd_scatterm_all                Calculates correction for downwelling scattered atmospheric radiation [MW 2012], section V.
!
! 6.  fd_tcos_eff                    Calculates effective cold space temperature taking into account the deviation between Rayleigh-Jeans and Planck law
!                                    as function of frequency [MWD 2012], Appendix D.
!
! 7.  get_emiss_wind                 Calculates isotropic wind induced emissivity at tabulated frequency values and reference EIA of 55.2. 
!                                    [MW 2012], section IV, Table 2.
!                               
! 8.  get_aharm_phir                 Calulates harmonic ocefficients of wind direction signal at tabulated frequency values and reference EIA of 55.2.
!                                    [MW 1012], section VI, Tables 3 + 4.
!
! 9.  get_aharm_phir_nad             Calculates harmonic ocefficients of wind direction signal at tabulated frequency values at nadir [MW 2012, section VI, eq. (26)]
!
! 10. get_sst_fac                    Precomputed value of E0(SST)/E0(SST_ref=20C) for faster computation [MW 2012], eq. (15) 
!
! 11. fd_xmea_win                    Provides wind speed polynomials for [MW 2012], eqs. (14) + (25) 
!

module RSS_RTM
save

public  ::    find_surface_tb, dielectric_meissner_wentz, fdem0_meissner_wentz, fd_emiss, fd_scatterm_all, fd_tcos_eff  
private ::    get_emiss_wind, get_aharm_phir, get_aharm_phir_nad, get_sst_fac, fd_xmea_win 

! external files 
! in ASCII     

! [MW 2012] Table 2
character(len=200), parameter     :: file_coeffs_wind_isotropic_A = 'finetune_emiss_wind.txt' 

! [MW 2012] Tables 3 + 4
character(len=200), parameter     :: file_coeffs_wind_direction_A = 'fit_emiss_phir_wind.txt' 

! [MW 2012] Tables 8, 9, 10, 11, 12
character(len=200), parameter     :: file_coeffs_sctterm_A        = 'mk_scatterm_table_all.txt' 

! [MW 2012] eq. (15) for fast computation
character(len=200), parameter     :: file_em0_ref_freq_sst_A      = 'em0_ref_freq_sst.txt' 


contains

! Master interface routine
subroutine find_surface_tb ( freq,tht,surtep,sal,ssws,phir,tran,tbdw,tc,        &
                             e0,ewind,omega,edirstokes,eharm,tbscat,tbsurf)

!   input 
!   name   descritpion                                            type      dimension     unit         established physical range                                       Default value              

!   freq   frequency                                              real(4)   scalar        GHz          [1.4, 89.0] for dielectric model + spec. emissivity    
!                                                                                                      [6.5, 89.0] for wind induced emissivity component
!                                                                                                      to be updated for L-band     

!   tht    earth incidence angle                                  real(4)   scalar        deg          [0,65]

!   surtep sea surface temperature                                real(4)   scalar        Kelvin       [271.15,313.15]                                   

!   sal    sea surface salinity                                   real(4)   scalar        ppt          [0,40]                                               optional    0.0 

!   ssws   sea surface wind speed                                 real(4)   scalar        m/s          [0,40]                                               optional    35.0                                           

!   phir   sea surface wind direction relative to azimuthal look  real(4)   scalar        deg          [0,360[                                              optional
!   upwind=0, downwind=180 deg

!   tran   atmospheric trnasmittance                              real(4)   scalar                     [0,1]                                                optional

!   tbdw   downwelling atmospheric brightness temperature         real(4)   scalar        Kelvin       >=0                                                  optional

!   tc     cold space temperature                                 real(4)   scalar        Kelvin       >=0                                                  optional   value from fd_tcos_eff  
!
!   output
!   name    descritpion                                           type      dimension     unit         physical range              

!   e0      specular sea surface emissivity  (1)=v, (2)=h         real(4)   vector(2)                  [0,1]                                                optional 

!   ewind   isotropic wind induced sea surface emissivity         real(4)   vector(2)                  [0,1]                                                optional
!           (1)=v, (2)=h


!   omega   omega term (correction for scattered downwelling      real(4)   vector(2)                  [0,1]                                                optional
!           reflection to include reflection from 
!           non specular directions)
!           calculated from geometric optics model using 
!           updated slope distribution
!           see amsr atbd  



!   edirstokes wind directional emissivity signal                 real(4)   vector(4)                  [0,1]                                                optional
!           modified stokes vector
!           (1)=v, (2)=h, (3)=s3, (4)=s4         


!   eharm harmonic coefficients of wind direction signal          real(4)   array(2,4)                 [0,1] 
!           modified stokes vector
!           index 1: 1=1st 2=2nd harmonic                          
!           index 2: (1)=v, (2)=h, (3)=s3, (4)=s4         
             

!   tbscat  scattered downweling radiation                        real(4)   vector(4)     Kelvin       >=0                                                  optional
!           modified stokes vector
!           (1)=v, (2)=h, (3)=s3, (4)=s4  

!   tbsurf  total brightness temperature from sea surface         real(4)   vector(4)     Kelvin       >=0                                                  optional        
!           emitted and reflected
!           modified stokes vector
!           (1)=v, (2)=h, (3)=s3, (4)=s4  


!   the following polarization basis vector convention for modified stokes parameters is used:
!   h = (k cross n) / abs(k cross n) (k: propagation direction from E to S/S, 
!   n: local mean sea surface normal upward from earth)                                        
!   v = h cross k
!   p = (v +  h)/sqrt(2)
!   m = (v -  h)/sqrt(2)
!   r = (v + ih)/sqrt(2) (ieee convention)
!   l = (v - ih)/sqrt(2) (ieee convention)
!   s3 = p - m  
!   s4 = l - r

implicit none
save

    real(4), intent(in)                               ::    freq
    real(4), intent(in)                               ::    tht
    real(4), intent(in)                               ::    surtep
    real(4), intent(in), optional                     ::    sal
    real(4), intent(in), optional                     ::    ssws
    real(4), intent(in), optional                     ::    phir
    real(4), intent(in), optional                     ::    tran
    real(4), intent(in), optional                     ::    tbdw
    real(4), intent(in), optional                     ::    tc

    real(4), intent(out), dimension(2),   optional	  ::    e0 
    real(4), intent(out), dimension(2),   optional	  ::    ewind 
    real(4), intent(out), dimension(2),   optional	  ::    omega 
    real(4), intent(out), dimension(4),   optional	  ::    edirstokes
    real(4), intent(out), dimension(2,4), optional	  ::    eharm
    real(4), intent(out), dimension(4),   optional	  ::    tbscat 
    real(4), intent(out), dimension(4),   optional	  ::    tbsurf 


    real(4)			                                  ::    xtc, xomegabar, sst, xsal, xssws, xphir
    real(4), dimension(2)                             ::    xe0, xewind, xscat, xomega
    real(4)				                              ::    costht, path, opacty
    real(4), dimension(4)                             ::    xestokes, xtbscat, xetot, xrtot
    real(4), dimension(2,4)                           ::    xeharm
    
    
    sst=surtep-273.15

    if (present(ssws)) then
      xssws=ssws
    else
      xssws=0.0
    endif

    if (present(phir)) then
      xphir=phir
    else
      xphir=-999.0
    endif


    if (present(sal)) then
      xsal=sal
    else
      xsal=35.0
    endif 

    ! effective cold space temperature including deviation from Rayleigh-Jeans law
    if (present(tc)) then
      xtc=tc
    else
      call fd_tcos_eff(freq, xtc)
    endif


    call fd_emiss(freq=freq,tht=tht,sst=sst,wind=xssws,phir=xphir, &
                  emiss_0=xe0, emiss_wind=xewind, emiss_phi=xestokes, eharm=xeharm)


    if (present(e0))           e0=xe0
    if (present(ewind))        ewind=xewind
    if (present(edirstokes))   edirstokes=xestokes
    
    if (present(eharm))        eharm=xeharm

    if (present(tbscat) .or. present(tbsurf) .or. present(omega)) then

        if ( .not.(present(tran)))   stop  ' need tran for computing tbscat. pgm stopped.'
        if ( .not.(present(tbdw)))   stop  ' need tbdw for computing tbscat. pgm stopped.'

        costht=cosd(tht)
        path=1.00035/sqrt(costht*costht+7.001225e-4)   
        ! (1+hratio)/sqrt(costht**2+hratio*(2+hratio)), hratio=.00035 [MW 2012] eq. (3)
        
        opacty=-alog(tran)/path
        call fd_scatterm_all(freq,tht,xssws,opacty, xscat)
        xomega(1:2) = xscat(1:2) /(tbdw + tran*xtc - xtc) ! [MW 2012] eq. (21)

        xetot(1:4) = (/xe0(1),xe0(2),0.0,0.0/) + (/xewind(1),xewind(2),0.0,0.0/) + xestokes(1:4)
        ! limit emissivity to [0,1]
        where(xetot(1:2)>=1.0) xetot(1:2)=1.0
        where(xetot(1:2)<=0.0) xetot(1:2)=0.0
      
        xrtot(1:2) = 1.0-xetot(1:2)
        xrtot(3:4) =    -xetot(3:4)
        xomegabar= (xomega(1)*xrtot(1) + xomega(2)*xrtot(2))/(xrtot(1)+xrtot(2))
        xtbscat(1:2) = ((1.0+xomega(1:2))*(tbdw+tran*xtc) - xomega(1:2)*xtc) * xrtot(1:2) 
        xtbscat(3:4) = ((1.0+xomegabar ) *(tbdw+tran*xtc) - xomegabar  *xtc) * xrtot(3:4)   
      
        if (present(tbscat)) tbscat = xtbscat
        if (present(tbsurf)) tbsurf = xetot*surtep + xtbscat
        if (present(omega))  omega=xomega

    endif
    
return
end subroutine find_surface_tb 


 
subroutine fd_emiss(freq,tht,sst,wind,phir,sal,     emiss_0,emiss_wind,emiss_phi,emiss_tot,eharm)
!  Calculates emissivity of specular surface (v/h), isotropic wind induced emissivity (v/h) and wind direction signal (v/h/S3/S4) 
!  [MW 2012, sections IV + VI.]
implicit none

     real(4), intent(in)                    :: freq,tht,sst
     real(4), intent(in), optional          :: wind
     real(4), optional, intent(in)          :: sal
     real(4), optional, intent(in)          :: phir

     
     real(4), optional, intent(out)         :: emiss_0(2),emiss_tot(2),emiss_wind(2),emiss_phi(4),eharm(2,4)

     real(4)                                :: xemiss_tot(2),xemiss(2),xemiss_phi(4),xsal 

     integer(4), parameter                  :: nstoke=4


     integer(4)                             :: ifreq1,ifreq2,ipol,istoke,iharm
     real(4)                                :: em0(2)
     real(4)                                :: qtht,wt,emiss1(2),emiss2(2),enad,h1,h2
     real(4)                                :: xphir
     real(4), save                          :: cos1phi,cos2phi,sin1phi,sin2phi
     real(4)                                :: aharm1(2,nstoke),aharm2(2,nstoke),aharm(2,nstoke),amp1,amp2,amp,anad(2,nstoke)

     real(4), save                          :: phirsv = 1.e30
     real(4), parameter                     :: thtref = 55.2	 
     real(4), dimension(2), parameter       :: xexp =(/4.,1.5/) ![MW 2012] section IV C.
     real(4), parameter, dimension(2,nstoke):: xexp_phir=reshape((/ 2.,2., 1.,4., 1.,4., 2.,2./), (/2,4/)) ![MW 2012] Table 5
     real(4), parameter, dimension(6)       :: freq0 =(/ 6.8,  10.7,  18.7,  23.8,  37.0, 85.5/)  
     !now referenced to windsat and ssmi

     if (present(phir)) then
        xphir=phir
     else
        xphir=-999.0
     endif

     if (present(sal)) then
        xsal=sal
     else
        xsal=35.0
     endif  

     if(freq.lt.6.5) stop 'pgm stopped, freq too small in oob in	fd_wind_emiss'

     qtht=tht
	 if(qtht.gt.65) qtht=65.  !qtht is just used for extrpolation for tht>thtref and i limit it to 60 deg

     if (present(emiss_0) .or. present(emiss_tot)) then
        call fdem0_meissner_wentz(freq=freq,tht=tht,sst=sst,salinity=xsal, em0=em0) 
     endif
      

    ifreq1=1
    if(freq.gt.freq0(2)) ifreq1=2
    if(freq.gt.freq0(3)) ifreq1=3  !between 18.7 and 37 ghz
    if(freq.gt.freq0(5)) ifreq1=5
 
    if(ifreq1.ne.3) then
        ifreq2=ifreq1+1
    else
        ifreq2=ifreq1+2
    endif
 
    wt=(freq-freq0(ifreq1))/(freq0(ifreq2)-freq0(ifreq1))
    if(freq.gt.freq0(ifreq2)) wt=1  !only occurs for freq>85.5

!   isotropic wind-induced emissivity
    if (present(emiss_wind) .or. present(emiss_tot)) then
        if(.not. (present(wind))) stop ' need wind speed for computing wind induced emissivity. pgm stopped.'
        call get_emiss_wind(ifreq1,sst,wind, emiss1)
        call get_emiss_wind(ifreq2,sst,wind, emiss2)
        xemiss=(1-wt)*emiss1 + wt*emiss2 
        !emiss is thtref value interpolated to input freq

        enad=0.5*(xemiss(1)+xemiss(2))

        do ipol=1,2
            if(tht.le.thtref) then
                xemiss(ipol)=enad         + (xemiss(ipol)-enad)*( tht/thtref)**xexp(ipol)
            else
                xemiss(ipol)=xemiss(ipol) + (xemiss(ipol)-enad)*(qtht-thtref)*xexp(ipol)/thtref  
            endif
        enddo  !ipol

    endif ! present isotropic wind induced emissivity

!   find emiss_phi
    if (present(emiss_phi) .or. present(emiss_tot) .or. present(eharm)) then
        if(.not.(present(wind))) stop ' need wind speed for computing wind induced emissivity. pgm stopped.'

        if(.not.present(eharm) .and.(xphir.lt.-998. .or. wind.le.3) ) then !-999. default for doing no correction
            xemiss_phi=0.0
        else  !find emiss_phi
            call get_aharm_phir(ifreq1,sst,wind, aharm1)  !aharm in terms of true stokes
            call get_aharm_phir(ifreq2,sst,wind, aharm2)  !aharm in terms of true stokes
            aharm=(1-wt)*aharm1 + wt*aharm2  !aharm is thtref value interpolated to input freq
            !     get nadir harmonic
            call get_aharm_phir_nad(ifreq1,freq0(ifreq1),sst,wind, amp1) 
            call get_aharm_phir_nad(ifreq2,freq0(ifreq2),sst,wind, amp2)  
            amp=(1-wt)*amp1 + wt*amp2
            anad=0  !most elements are zero
            anad(2,2)=  amp
            anad(2,3)= -amp
            do istoke=1,nstoke
            do iharm=1,2
                if(tht.le.thtref) then
                    aharm(iharm,istoke)= anad(iharm,istoke) +(aharm(iharm,istoke)-anad(iharm,istoke))* (tht/thtref)**xexp_phir(iharm,istoke)
                else
                    aharm(iharm,istoke)=aharm(iharm,istoke) +(aharm(iharm,istoke)-anad(iharm,istoke))*(qtht-thtref)*xexp_phir(iharm,istoke)/thtref  
                endif
            enddo  !iharm
            enddo  !istoke
            !     convert back from true stokes to v and h
            do iharm=1,2
                h1=aharm(iharm,1) + 0.5*aharm(iharm,2) !(v+h)/2 + (v-h)/2=v
                h2=aharm(iharm,1) - 0.5*aharm(iharm,2) !(v+h)/2 - (v-h)/2=h
                aharm(iharm,1)=h1
                aharm(iharm,2)=h2
            enddo

            if(abs(xphir-phirsv).gt.0.01) then
                phirsv=xphir
                cos1phi=cosd(  xphir)
                cos2phi=cosd(2*xphir)
                sin1phi=sind(  xphir)
                sin2phi=sind(2*xphir)
            endif

            xemiss_phi(1:2)=aharm(1,1:2)*cos1phi + aharm(2,1:2)*cos2phi
            xemiss_phi(3:4)=aharm(1,3:4)*sin1phi + aharm(2,3:4)*sin2phi
            xemiss_phi(4) = - xemiss_phi(4) ! IEEE convention

        endif ! (xphir.lt.-998. .or. wind.le.3)
        
    endif  ! emiss_phi

    if (present(emiss_tot))  then
        xemiss_tot = em0 + xemiss  + xemiss_phi(1:2)
        emiss_tot  = xemiss_tot
    endif
 
    if (present(emiss_0))    emiss_0   =em0
    if (present(emiss_wind)) emiss_wind=xemiss
    if (present(emiss_phi))  emiss_phi =xemiss_phi

    if (present(eharm)) then
        eharm=aharm
        eharm(:,4)=-aharm(:,4) ! IEEE
        if (wind<=3.0) eharm=0.0
    endif

return
end subroutine fd_emiss



subroutine get_emiss_wind(ifreq,sst,wind,   emiss)
! [MW 2012] section IV
implicit none

    integer(4), intent(in)                ::  ifreq
    real(4), intent(in)                   ::  sst,wind
    real(4), dimension(2), intent(out)    ::  emiss

    real(4), save                         ::  acoef(5,2,6)
    real(4)                               ::  sst_fac(2)
    real(8)                               ::  xmea(5)
    
    integer(4)                            ::  i,j,k
    integer(4)                            ::  ii,jj,kk


    integer(4), save :: istart=1

    if(istart.eq.1) then
      istart=0
      
      !open(unit=3,file=file_coeffs_wind_isotropic,form='binary',status='old',action='read')
      !read(3) acoef
      !close(3)
    
      open(unit=3,file=file_coeffs_wind_isotropic_A,form='formatted',action='read',status='old')
      do i=1,5
      do j=1,2
      do k=1,6
      
      read(3,6003) ii,jj,kk,acoef(i,j,k)
      6003 format(1x,i3,1x,i3,1x,i3,1x,e20.8)
      
      enddo
      enddo
      enddo  
    
      close(3)
    
    
    endif

    if(ifreq.eq.4) stop 'ifreq oob in get_emiss_wind, pgm stopped'

    call  fd_xmea_win(wind, xmea)
    emiss(1)=dot_product(acoef(:,1,ifreq),xmea) 
    emiss(2)=dot_product(acoef(:,2,ifreq),xmea) 

    call get_sst_fac(ifreq,1,sst, sst_fac(1))
    call get_sst_fac(ifreq,2,sst, sst_fac(2))

    emiss=emiss*sst_fac


return
end subroutine get_emiss_wind




subroutine get_aharm_phir(ifreq,sst,wind,   aharm)
! [MW 2012] section VI
implicit none

    integer(4), parameter                           :: nstoke=4
    
    integer(4), intent(in)                          :: ifreq 
    real(4), intent(in)                             :: sst, wind
    
    real(4), dimension(2,nstoke), intent(out)       :: aharm

    integer(4)                                      :: istoke,iharm
    
    real(4)                                         :: h1,h2
    real(4), save                                   :: bcoef(5,2,nstoke,6)
    real(4)                                         :: sst_fac(nstoke)
    real(8)                                         :: xmea(5)

    integer(4)                                      :: i,j,k,l
    integer(4)                                      :: ii,jj,kk,ll


    integer(4), save  :: istart=1
 
    if(istart.eq.1) then
        istart=0
        
        !open(unit=3,file=file_coeffs_wind_direction,form='binary',status='old',action='read')
        !read(3) bcoef
        !close(3) 
        
        open(unit=3,file=file_coeffs_wind_direction_A,form='formatted',action='read',status='old')
        do i=1,5
        do j=1,2
        do k=1,nstoke
        do l=1,6
        
        read(3,6004) ii,jj,kk,ll,bcoef(i,j,k,l)
        6004 format(1x,i3,1x,i3,1x,i3,1x,i3,1x,e20.8)
        
        enddo
        enddo
        enddo
        enddo        
        
        close(3)
        
    endif

    if(ifreq.eq.4) stop 'ifreq oob in get_aharm_phir, pgm stopped'
    
    call  fd_xmea_win(wind, xmea)

    call get_sst_fac(ifreq,1,sst,   sst_fac(1))
    call get_sst_fac(ifreq,2,sst,   sst_fac(2))
    sst_fac(3:4)=0.5*(sst_fac(1)+sst_fac(2))

    do istoke=1,nstoke
        aharm(1,istoke)=sst_fac(istoke)*dot_product(xmea,bcoef(:,1,istoke,ifreq))
        aharm(2,istoke)=sst_fac(istoke)*dot_product(xmea,bcoef(:,2,istoke,ifreq))
    enddo  !istoke

!     convert to true stokes paramters,ie (v+h)/2 and v-h, rather than v and h in order to do tht adjustment
    do iharm=1,2
        h1=0.5*(aharm(iharm,1)+aharm(iharm,2))
        h2=     aharm(iharm,1)-aharm(iharm,2)
        aharm(iharm,1)=h1
        aharm(iharm,2)=h2
    enddo
  
return
end subroutine get_aharm_phir



subroutine get_aharm_phir_nad(ifreq,freq,sst,wind,      amp)
! [MW 2012], section VI] eq. (26).
implicit none

    integer(4), intent(in)    :: ifreq
    real(4), intent(in)       :: freq,sst,wind
    real(4), intent(out)      :: amp
    
    real(4)                   :: amp_10_nad,ywind,qfreq
    real(4)                   :: sst_fac

    qfreq=freq
    if(qfreq.gt.37) qfreq=37

    if(freq.lt.3) then
        amp_10_nad=.2/290.
    else
        amp_10_nad=2*(1. - 0.9*alog10(30./qfreq))/290.
    endif
 
    ywind=wind
    if(wind.lt. 0) ywind= 0
    if(wind.gt.15) ywind=15

    amp=amp_10_nad*ywind*(ywind - ywind**2/22.5)/55.5556
    call get_sst_fac(ifreq,0,sst, sst_fac)
    amp=amp*sst_fac

return
end subroutine get_aharm_phir_nad



subroutine get_sst_fac(ifreq,ipol,sst,      sst_fac)
! [MW 2012], eq. (15)
!ipol=0 denotes nadir value, sst_fac=em0(sst)/em0(sst=20)
!rcoef values for the ratio of nadir em0(sst)/em0(sst=20) were computed offline 
!so that they are available for fast computation

implicit none

    integer(4), intent(in)  :: ifreq,ipol
    real(4), intent(in)     :: sst
    
    real(4), intent(out)    :: sst_fac

    integer(4), save        :: istart=1
    real(4)                 :: xmea(3)
    real(4), save           :: rcoef(3,0:2,6)

    integer(4)              :: i,j,k
    integer(4)              :: ii,jj,kk


    if(istart.eq.1) then
        istart=0
        
        !open(unit=3,file=file_em0_ref_freq_sst,status='old',form='binary',action='read')
        !read(3) rcoef
        !close(3)
    
        open(unit=3,file=file_em0_ref_freq_sst_A,form='formatted',action='read',status='old')
        
        do i=1,3 
        do j=0,2
        do k=1,6
        
        read(3,6001) ii,jj,kk,rcoef(i,j,k)
        6001 format(1x,i3,1x,i3,1x,i3,1x,e20.8)
        
        enddo
        enddo
        enddo
        close(3)
    
    
    endif

    xmea(1)= sst-20
    xmea(2)=xmea(1)*xmea(1)
    xmea(3)=xmea(1)*xmea(2)
    sst_fac=1 + dot_product(rcoef(:,ipol,ifreq),xmea)

return
end subroutine get_sst_fac



subroutine fd_xmea_win(wind,        xmea)
! Provides wind speed polynomials for [MW 2012, eqs. (14) + (25)] 
implicit none
 
    real(4), intent(in)                     :: wind
    real(8), dimension(5), intent(out)      :: xmea
    
    real(4)                                 :: x,dif
    real(4), parameter                      :: wcut =20.
    
    x=wind
    if(x.lt.0) x=0
    
    xmea(1)=x
    if(x.le.wcut) then
        xmea(2)=xmea(1)*x
        xmea(3)=xmea(2)*x
        xmea(4)=xmea(3)*x
        xmea(5)=xmea(4)*x
    else
        dif=x-wcut
        xmea(2)=2*dif*wcut       + wcut**2
        xmea(3)=3*dif*wcut**2    + wcut**3
        xmea(4)=4*dif*wcut**3    + wcut**4
        xmea(5)=5*dif*wcut**4    + wcut**5
    endif
 
return
end subroutine fd_xmea_win



subroutine fd_scatterm_all(freq,tht,wind,opacty,    xscat)
! [MW 2012], section V.
implicit none

    real(4), intent(in)                             :: freq,tht,wind,opacty
    real(4), dimension(2), intent(out)              :: xscat
    
    real(4)                                         :: xlog_freq,xscat1(2),xscat2(2)
    real(4)                                         :: a1,a2,b1,b2,c1,c2,brief,d1,d2
    real(4), save                                   :: scatterm(91,50,26,13,2)

    integer(4), save                                :: istart=1
    integer(4)                                      :: i1,i2,j1,j2,k1,k2,l1,l2

    integer(4)                                      :: i,j,k,l,m
    integer(4)                                      :: ii,jj,kk,ll,mm


    if(istart.eq.1) then
        
        istart=0
        
        !open(unit=3,file=file_coeffs_sctterm,status='old',form='binary',action='read')
        !read(3) scatterm
        !close(3)
        
        open(unit=3,file=file_coeffs_sctterm_A,form='formatted',action='read',status='old')
        do i=1,91
        do j=1,50
        do k=1,26
        do l=1,13
        do m=1,2
        
        read(3,6002) ii,jj,kk,ll,mm,scatterm(i,j,k,l,m)       
        6002 format(1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,e20.8)
        
        enddo
        enddo
        enddo
        enddo
        enddo
         
        close(3)
        
    endif

    ! check inputs

    if(freq.lt.1 .or. freq.gt.200) stop 'freq oob in fd_scatterm, pgm stopped'
    if(tht .lt.0 .or.  tht.gt. 90) stop 'tht  oob in fd_scatterm, pgm stopped'
    if(wind.lt.0 .or. wind.gt.100) stop 'wind oob in fd_scatterm, pgm stopped'
    if(opacty.lt.0)                stop 'opacty oob in fd_scatterm, pgm stopped'

    xlog_freq=alog10(freq)
    
    ! multi-linear interpolation from table values
    
    brief=tht
    if(brief.gt.89.99) brief=89.99
    i1=1+brief
    i2=i1+1
    a1=i1-brief
    a2=1.-a1
    
    brief=wind
    if(brief.gt.24.99) brief=24.99
    j1=1+brief
    j2=j1+1
    b1=j1-brief
    b2=1-b1
    
    brief=xlog_freq/0.2
    if(brief.gt.11.99) brief=11.99
    k1=1+brief
    k2=k1+1
    c1=k1-brief
    c2=1-c1
    
    brief=opacty/0.025
    if(brief.gt.48.99) brief=48.99
    l1=1+brief
    l2=l1+1
    d1=l1-brief
    d2=1-d1
    
    xscat1= &
    a1*b1*(c1*scatterm(i1,l1,j1,k1,:)+c2*scatterm(i1,l1,j1,k2,:))+ &
    a1*b2*(c1*scatterm(i1,l1,j2,k1,:)+c2*scatterm(i1,l1,j2,k2,:))+ &
    a2*b1*(c1*scatterm(i2,l1,j1,k1,:)+c2*scatterm(i2,l1,j1,k2,:))+ &
    a2*b2*(c1*scatterm(i2,l1,j2,k1,:)+c2*scatterm(i2,l1,j2,k2,:))

    xscat2= &
    a1*b1*(c1*scatterm(i1,l2,j1,k1,:)+c2*scatterm(i1,l2,j1,k2,:))+ &
    a1*b2*(c1*scatterm(i1,l2,j2,k1,:)+c2*scatterm(i1,l2,j2,k2,:))+ &
    a2*b1*(c1*scatterm(i2,l2,j1,k1,:)+c2*scatterm(i2,l2,j1,k2,:))+ &
    a2*b2*(c1*scatterm(i2,l2,j2,k1,:)+c2*scatterm(i2,l2,j2,k2,:))
    
    xscat=d1*xscat1 + d2*xscat2
 
return
end subroutine fd_scatterm_all


subroutine fdem0_meissner_wentz(freq,tht,sst,salinity,      em0)
! Compute specular emissivity using Frsenel equations and MW dielectric model.
implicit none

    real(4), intent(in)                         :: freq,tht,sst,salinity

    real(4), dimension(2), intent(out)          :: em0
    
    real(4), parameter                          :: f0=17.97510
 

    real(4)                                     :: costht,sinsqtht
    real(4)                                     :: e0s,e1s,e2s,n1s,n2s,sig
    
    complex(4)                                  :: permit,esqrt,rh,rv
    complex(4), parameter                       :: j=(0.,1.)
    
 
    call dielectric_meissner_wentz(sst,salinity,  e0s,e1s,e2s,n1s,n2s,sig)

    costht=cosd(tht)
    sinsqtht=1.-costht*costht


!   debye law (2 relaxation wavelengths)
    permit = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s + j*sig*f0/freq
    permit = conjg(permit)
    
    esqrt=csqrt(permit-sinsqtht)
    rh=(costht-esqrt)/(costht+esqrt)
    rv=(permit*costht-esqrt)/(permit*costht+esqrt)
    em0(1)  =1.-rv*conjg(rv)
    em0(2)  =1.-rh*conjg(rh)
 
return
end subroutine fdem0_meissner_wentz

 


subroutine dielectric_meissner_wentz(sst_in,s,   e0s,e1s,e2s,n1s,n2s,sig)
!
!     complex dielectric constant: eps
!     [MW 2004, MW 2012, MWR 2014].
!     
!     Changes from [MW 2012]:
!     1. Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
!     2. Changed SST behavior of coefficient b2 from:
!     b2 = 1.0 + s*(z(10) + z(11)*sst) to
!     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
!
!!
!     input:
!     name   parameter  unit  range
!     sst      sst        [c]   -25 c to 40 c for pure water
!                               -2  c to 34 c for saline water
!     s      salinity   [ppt]  0 to 40
!
!     output:
!     eps    complex dielectric constant
!            negative imaginary part to be consistent with wentz1 convention
!

implicit none


    real(4), intent(in)  :: sst_in,s
    real(4), intent(out) :: e0s,e1s,e2s,n1s,n2s,sig
 
    real(4), dimension(11), parameter :: &
      x=(/ 5.7230e+00, 2.2379e-02, -7.1237e-04, 5.0478e+00, -7.0315e-02, 6.0059e-04, 3.6143e+00, &
           2.8841e-02, 1.3652e-01,  1.4825e-03, 2.4166e-04 /)
    
    real(4), dimension(13), parameter :: &
      z=(/ -3.56417e-03,  4.74868e-06,  1.15574e-05,  2.39357e-03, -3.13530e-05, &
            2.52477e-07, -6.28908e-03,  1.76032e-04, -9.22144e-05, -1.99723e-02, &
            1.81176e-04, -2.04265e-03,  1.57883e-04  /)  ! 2004

    real(4), dimension(3), parameter :: a0coef=(/ -0.33330E-02,  4.74868e-06,  0.0e+00/)
    real(4), dimension(5), parameter :: b1coef=(/0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08/)
 
    real(4) :: e0,e1,e2,n1,n2
    real(4) :: a0,a1,a2,b1,b2
    real(4) :: sig35,r15,rtr15,alpha0,alpha1

    real(4) :: sst,sst2,sst3,sst4,s2
    
    sst=sst_in
    if(sst.lt.-30.16) sst=-30.16  !protects against n1 and n2 going zero for very cold water
    
    sst2=sst*sst
    sst3=sst2*sst
    sst4=sst3*sst

    s2=s*s
 
    !     pure water
    e0    = (3.70886e4 - 8.2168e1*sst)/(4.21854e2 + sst) ! stogryn et al.
    e1    = x(1) + x(2)*sst + x(3)*sst2
    n1    = (45.00 + sst)/(x(4) + x(5)*sst + x(6)*sst2)
    e2    = x(7) + x(8)*sst
    n2    = (45.00 + sst)/(x(9) + x(10)*sst + x(11)*sst2)
    
    !     saline water
    !     conductivity [s/m] taken from stogryn et al.
    sig35 = 2.903602 + 8.60700e-2*sst + 4.738817e-4*sst2 - 2.9910e-6*sst3 + 4.3047e-9*sst4
    r15   = s*(37.5109+5.45216*s+1.4409e-2*s2)/(1004.75+182.283*s+s2)
    alpha0 = (6.9431+3.2841*s-9.9486e-2*s2)/(84.850+69.024*s+s2)
    alpha1 = 49.843 - 0.2276*s + 0.198e-2*s2
    rtr15 = 1.0 + (sst-15.0)*alpha0/(alpha1+sst)
    
    sig = sig35*r15*rtr15
    
    !    permittivity
    a0 = exp(a0coef(1)*s + a0coef(2)*s2 + a0coef(3)*s*sst)  
    e0s = a0*e0
    
    if(sst.le.30) then
        b1 = 1.0 + s*(b1coef(1) + b1coef(2)*sst + b1coef(3)*sst2 + b1coef(4)*sst3 + b1coef(5)*sst4)
    else
        b1 = 1.0 + s*(9.1873715e-04 + 1.5012396e-04*(sst-30))
    endif
      
    n1s = n1*b1
    
    a1  = exp(z(7)*s + z(8)*s2 + z(9)*s*sst)
    e1s = e1*a1

    b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30))
    n2s = n2*b2
    
    a2 = 1.0  + s*(z(12) + z(13)*sst)
    e2s = e2*a2
    
return
end subroutine  dielectric_meissner_wentz

 

subroutine fd_tcos_eff(freq,    tcos_eff)
!     Calculates effective cold space temperature taking into account 
!     the deviation between Rayleigh-Jeans and Planck law
!     as function of frequency [MWD 2012], Appendix D.
!   
!     for these routine the term b is the flux 
!    (2*h*f**3/(c**2*(dexp(h*f/(k*t))-1))) divided by  2*k*f**2/c**2
implicit none


    real(4), intent(in)         :: freq
    real(4), intent(out)        :: tcos_eff


    real(8), parameter          :: tcos=2.73
    real(8), parameter          :: teff=63.  !selected to provide optimum fit over 60-300 k range
    real(8), parameter          :: h=6.6260755d-34
    real(8), parameter          :: k= 1.380658d-23
    real(8), parameter          :: a=h/k
    real(8)                     :: x,b1,b2

    x=a*freq*1.d9
    b1=x/(dexp(x/tcos)-1)
    b2=x/(dexp(x/teff)-1)
    tcos_eff=b1-b2+teff

return
end subroutine fd_tcos_eff



end module RSS_RTM