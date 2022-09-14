include 'RSS_RTM_module.f90'  
!
!  RSS_RTM_driver.f90 
!
!
!
!****************************************************************************
!
!  PROGRAM: RSS_RTM_driver
!
!  PURPOSE:  driver routine for RSS RTM
!  this routine creates the test output file "control_output.txt"
!
!  Reference 
!  [MW 2012]:
!  T. Meissner, and F. Wentz, 
!  The Emissivity of the Ocean Surface between 6 - 90 GHz over a Large Range of Wind Speeds and Earth Incidence Angles, 
!  IEEE Transactions on Geoscience and Remote Sensing, vol 50 (8), 3004-3026, doi: 10.1109/TGRS.2011.2179662, 2012. 
!
!****************************************************************************

program RSS_RTM_driver
use RSS_RTM
implicit none

    integer(4)                       :: i
    real(4)                          :: freq, tht, sst, sal, wind, surtep, tran,  td, tbdw
    real(4), dimension(2)            :: ewind   
    real(4), dimension(2)            :: omega    
    real(4), dimension(2,4)          :: eharm    
    real(4), parameter               :: teff = 290.0
    character(len=200), parameter    :: control_output='control_output.txt'
   
    open(unit=4,file=control_output,action='write',form='formatted')

    write(4,*) 'Contol Output of RSS RTM module'    
    write(4,*)
    
    
    ! Case 1: MW 2012 Figure 8
    freq=37.0
    tht=55.2
    sal=35.0
    sst=20.0
    surtep=sst+273.15
    write(4,*)
    write(4,*) ' [MW 2012] Figure 8'
    write(4,*) ' FREQUENCY = ', freq,' GHz'    
    write(4,*) ' EIA = ', tht,' deg'    
    do i=0,25
        wind=float(i)
        call find_surface_tb (freq=freq,tht=tht,surtep=surtep,sal=sal,ssws=wind,    ewind=ewind)
        write(4,101) wind,  ewind(1:2)*teff
        101 format(1x,f4.1,2x,f10.3,1x,f10.3)
    enddo
    
    
    ! Case 2: MW 2012 Figure 5
    freq=18.7
    sal=35.0
    sst=20.0
    surtep=sst+273.15
    wind=12.1
    write(4,*)
    write(4,*) ' [MW 2012] Figure 5'
    write(4,*) ' FREQUENCY = ', freq,' GHz'    
    write(4,*) ' WIND= ', wind,' m/s'    
    do i=0,6
        tht=float(i)*10.0
        call find_surface_tb (freq=freq,tht=tht,surtep=surtep,sal=sal,ssws=wind,    ewind=ewind)
        write(4,102) tht,  ewind(1:2)*teff
102     format(1x,f6.1,2x,f10.3,1x,f10.3)
    enddo
    

    ! Case 3: MW 2012 Figure 9
    tht=55.2
    sal=35.0
    sst=20.0
    surtep=sst+273.15
    wind=12.1
    write(4,*)
    write(4,*) ' [MW 2012] Figure 9'
    write(4,*) ' EIA = ', tht,' deg'    
    write(4,*) ' WIND= ', wind,' m/s'    
    do i=1,6
        
        select case (i)
        case (1) 
            freq=6.8
        case (2) 
            freq=10.7
        case (3) 
            freq=18.7
        case (4) 
            freq=23.8
        case (5) 
            freq=37.0
        case (6) 
            freq=85.5         
        end select

        call find_surface_tb (freq=freq,tht=tht,surtep=surtep,sal=sal,ssws=wind,    ewind=ewind)
        write(4,103) freq,  ewind(1:2)*teff
103     format(1x,f6.1,2x,f10.3,1x,f10.3)
    enddo
    
    ! Case 4: MW 2012 Figure 10
    freq=18.7
    tht=55.2
    sal=35.0
    wind=7.0
    sst=20.0
    surtep=sst+273.15
    td = 281.0
    write(4,*)
    write(4,*) ' [MW 2012] Figure 10'
    write(4,*) ' FREQUENCY = ', freq,' GHz'    
    write(4,*) ' EIA = ', tht,' deg'
    write(4,*) ' WIND ', wind,' m/s'
    do i=0,10
        tran=float(i)*0.1
        if (tran >= 0.9999999) tran=0.9999999
        tbdw = td*(1.0-tran)
        call find_surface_tb (freq=freq,tht=tht,surtep=surtep,sal=sal,ssws=wind,tran=tran,tbdw=tbdw,    omega=omega)
        write(4,104) tran, tbdw,omega(1:2)
104     format(1x,f4.1,1x,f8.3,5x,f8.3,1x,f8.3)
    enddo




    ! Case 5: MW 2012 Figures 12 + 13
    freq=10.7
    tht=50.1
    sal=35.0
    sst=20.0
    surtep=sst+273.15
    write(4,*)
    write(4,*) ' [MW 2012] Figures 12/13'
    write(4,*) ' FREQUENCY = ', freq,' GHz'    
    write(4,*) ' EIA = ', tht,' deg'    
    do i=0,25
        wind=float(i)
        call find_surface_tb (freq=freq,tht=tht,surtep=surtep,sal=sal,ssws=wind,    eharm=eharm)
        write(4,105) wind,  eharm(1:2,1)*teff, eharm(1:2,2)*teff, eharm(1:2,3)*teff, eharm(1:2,4)*teff
105     format(1x,f4.1,2x,4(2(f6.3,1x),1x))
    enddo


   close(4)


    
stop ' normal end'
end program RSS_RTM_driver

