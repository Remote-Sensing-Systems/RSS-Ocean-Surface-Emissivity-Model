# RSS-Ocean-Surface-Emissivity-Model
Remote Sensing Systems Ocean Surface Emissivity Model 6 - 90 GHz

 Remote Sensing System
 Santa Rosa, California
 
 Contact 
 Thomas Meissner
 meissner@remss.com

 September 14, 2022

 RSS_RTM_2012
 published on GitHub under MIT license 

 If the RTM plays an essential role in the generation of new research results, please acknowledge/credit RSS and let us know.



 FORTRAN 90 module for computing ocean surface emissivity at microwave frequencies
 Citation:
 T. Meissner and F.  Wentz, RSS Ocean Surface Emissivity Model, doi: 10.21982/M84S3Q.


 References:
 [MW 2004]:   T. Meissner and F. Wentz, 
              "The Complex Dielectric Constant of Pure and Sea Water from Microwave Satellite Observations", 
              IEEE Transactions on Geoscience and Remote Sensing, vol 42 (9), 1836-1849, doi: 10.1109/TGRS.2004.831888, 2004.               

 [MW 2012]:   T. Meissner, and F. Wentz, 
              "The Emissivity of the Ocean Surface between 6 - 90 GHz over a Large Range of 
              Wind Speeds and Earth Incidence Angles", 
              IEEE Transactions on Geoscience and Remote Sensing, vol 50 (8), 3004-3026, doi: 10.1109/TGRS.2011.2179662, 2012.     

 
 [MWR 2014]:  T. Meissner, F. Wentz, and L. Ricciardulli, 
              "The Emission and Scattering of L-band Microwave Radiation 
              from Rough Ocean Surfaces and Wind speed Measurements from Aquarius", 
              Journal of Geophysical Research: Oceans, 119, doi:10.1002/2014JC009837, 2014.    

 [MWD 2012]:  T. Meissner, F. Wentz and D. Draper, 
              "GMI Calibration Algorithm and Analysis Theoretical Basis Document", 
              report number 041912, Version-G, Remote Sensing Systems, Santa Rosa, CA, 124 pp.
              doi: 10.56236/RSS-au,    
              https://doi.org/10.56236/RSS-au 


 Content:
1. readme.mc:           This file
2. RSS_RTM__module.f90: FORTRAN 90 module for computing ocean surface emissivity at microwave frequencies.
3. RSS_RTM_driver.f90:  FORTRAN 90 sample driver.
4. control_output.txt:  Sample control output created by driver routine.
5. external_files:      4 external ASCII files (tables):
                        finetune_emiss_wind.txt
                        fit_emiss_phir_wind.txt   
                        mk_scatterm_table_all.txt
                        em0_ref_freq_sst.txt

6. meissner_wentz_TGRS_2012.pdf
[MW 2012]:   
T. Meissner, and F. Wentz, 
"The Emissivity of the Ocean Surface between 6 - 90 GHz over a Large Range of Wind Speeds and Earth Incidence Angles", 
IEEE Transactions on Geoscience and Remote Sensing, vol 50 (8), 3004-3026, doi: 10.1109/TGRS.2011.2179662, 2012. 


 Release Notes:

 The following changes/updates from [MW 2012] are contained in the code.
 1. dielectric_meissner_wentz: Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
 2. dielectric_meissner_wentz: Changed SST behavior of coefficient b2 from:
     b2 = 1.0 + s*(z(10) + z(11)*sst) to
     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
 3. Introduced SST dependence for wind direction signal similar than eq. (15) in [MW 2012].  


 Routines:

 1.  find_surface_tb                Master wrapper routine that caclulates all components of the ocean surface emissivity. Output user selected.

 2.  dielectric_meissner_wentz      Dielectric model of sea and pure water [MW 2004], with minor updates in [MW 2012].

 3.  fdem0_meissner_wentz           Caclulates emissivity of specular surface [MW 2004]

 4.  fd_emiss                       Calculates emissivity of specular surface (v/h), isotropic wind induced emissivity (v/h) and wind direction signal (v/h/S3/S4) 
                                    [MW 2012], sections IV + VI.

 5.  fd_scatterm_all                Calculates correction for downwelling scattered atmospheric radiation [MW 2012], section V.

 6.  fd_tcos_eff                    Calculates effective cold space temperature taking into account the deviation between Rayleigh-Jeans and Planck law
                                    as function of frequency [MWD 2012], Appendix D.

 7.  get_emiss_wind                 Calculates isotropic wind induced emissivity at tabulated frequency values and reference EIA of 55.2. 
                                    [MW 2012], section IV, Table 2.
                               
 8.  get_aharm_phir                 Calulates harmonic ocefficients of wind direction signal at tabulated frequency values and reference EIA of 55.2.
                                    [MW 1012], section VI, Tables 3 + 4.

 9.  get_aharm_phir_nad             Calculates harmonic ocefficients of wind direction signal at tabulated frequency values at nadir [MW 2012, section VI, eq. (26)]

 10. get_sst_fac                    Precomputed value of E0(SST)/E0(SST_ref=20C) for faster computation [MW 2012], eq. (15) 

 11. fd_xmea_win                    Provides wind speed polynomials for [MW 2012], eqs. (14) + (25) 
