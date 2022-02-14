% The Snow, Ice, and Aerosol Radiative (SNICAR) Model: stand-alone version
%
% Original version: written by Mark Flanner (flanner@umich.edu) Flanner et al. 2007 JGR
% Current version: Updated by Cenlin He (cenlinhe@ucar.edu) He et al. 2018 ACP
% Original version assumes aerosol externally mixed with spherical snow grains 
% Current version includes nonspherical snow grain & BC-snow internal mixing 
% based on parameterizations from [He et al. 2017,doi:10.1175/JCLI-D-17-0300.1]
%  
% This function calculates radiative fluxes at the interfaces of
% multiple snow layers using a 2-stream approximation with a
% tri-diagonal matrix solution, presented by Toon, O.B. et al.,
% Journal for Geophysical Research, vol 94 D13 p 16287, 1989.  Thus, a
% vertically inhomogenous medium can be represented.
%
% The function reads in Mie optical properties from external NetCDF
% files based on defined properties of each layer.
%
% Example input parameters are defined below. The user has option of
% running the program as a script. To do so, comment out the function
% line and set '1==1' below:
%
%
%%%%%%%%%%  Input parameters: %%%%%%%%%%%
% BND_TYP:      Spectral grid (=1 for 470 bands. This is the
%               only functional option in this distribution)
% DIRECT:       Direct or diffuse incident radiation (1=direct, 0=diffuse)
% APRX_TYP:     Two-Stream Approximation Type:
%                  1 = Eddington
%                  2 = Quadrature
%                  3 = Hemispheric Mean
%                  NOTE: Delta-Eddington Approximation is probably
%                  best for the visible spectrum, but can provide
%                  negative albedo in the near-IR under diffuse
%                  light conditions. Hence, the Hemispheric Mean
%                  approximation is recommended for general use.
% DELTA:        1=Use Delta approximation (Joseph, 1976), 0=don't use
% coszen:       cosine of solar zenith angle (only applies when DIRECT=1)
% R_sfc:        broadband albedo of underlying surface 
%                  (user can also define a spectrally-dependent albedo below)
% dz:           array of snow layer thicknesses [meters]. Top layer has index 1. 
%                  The length of this array defines number of snow layers
% rho_snw:      array of snow layer densities [kg m-3]. Must have same length as dz
% rds_snw:      array of snow layer effective grain radii [microns]. Must have same length as dz
% nbr_aer:      number of aerosol species in snowpack
% mss_cnc_sot1: mass mixing ratio of black carbon species 1 (uncoated BC)
%                 (units of parts per billion, ng g-1)
% mss_cnc_sot2: mass mixing ratio of black carbon species 2 (sulfate-coated BC)
%                 (units of parts per billion, ng g-1)
% mss_cnc_sot3: mass mixing ratio of black carbon species 3 (internally mixed with snow)
%                 (units of parts per billion, ng g-1) add by che
% mss_cnc_dst1: mass mixing ratio of dust species 1 (radii of 0.05-0.5um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst2: mass mixing ratio of dust species 2 (radii of 0.5-1.25um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst3: mass mixing ratio of dust species 3 (radii of 1.25-2.5um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_dst4: mass mixing ratio of dust species 4 (radii of 2.5-5.0um)
%                 (units of parts per billion, ng g-1)
% mss_cnc_ash1: mass mixing ratio of volcanic ash species 1
%                 (units of parts per billion, ng g-1)
% fl_sot1:      name of file containing optical properties for BC species 1
% fl_sot2:      name of file containing optical properties for BC species 2
% fl_dst1:      name of file containing optical properties for dust species 1
% fl_dst2:      name of file containing optical properties for dust species 2
% fl_dst3:      name of file containing optical properties for dust species 3
% fl_dst4:      name of file containing optical properties for dust species 4
% fl_ash1:      name of file containing optical properties for ash species 1
%
% NOTE:  The spectral distribution of incident radiation is, by
% default, one typical of a mid-latitude winter. The user can
% change the incident flux distribution by changing the text files
% loaded below under "Set incident solar flux spectral distribution:"
%
%
%%%%% che %%%%%% options added by Cenlin He, Apr 2018
% SNO_SHP:      Snow grain shape
%               1=sphere (SNICAR default shape)
%               2=spheroid; 3=hexagonal plate; 4=koch snowflake;
% SNO_fs:       Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
%               0=use recommended default value (He et al. 2017 JC);
%               others(0<fs<1)= use user-specified value
%               only activated when SNO_SHP > 1 (i.e. nonspherical)
% SNO_AR:       Aspect ratio: ratio of grain width to length
%               0=use recommended default value (He et al. 2017 JC);
%               others(0.1<fs<20)= use user-specified value
%               only activated when SNO_SHP > 1 (i.e. nonspherical)
% BC_SNO_MIX:   BC-snow mixing option
%               1= only BC-snow external mixing (SNICAR default)
%               2= BC-snow internal + external mixing
%%%%% che %%%%%%

%%%%%  Output data: %%%%%
% All output data is contained in the matrix "data_out"
% data_out(:,1) = wavelength (microns, m^-6)
% data_out(:,2) = spectrally-dependent albedo;
% data_out(1,3) = broadband (0.3-5.0um) solar albedo
% data_out(2,3) = visible (0.3-0.7um) albedo;
% data_out(3,3) = near-IR (0.7-5.0um) albedo;
% data_out(4,3) = total solar absorption by snowpack (W/m2)
%                 (not including underlying substrate)
% data_out(5,3) = solar absorption in top snow layer (W/m2)
% data_out(6,3) = solar absorption in second snow layer (W/m2)
% (see "data_out" fields near end of script for more output and to
% add new output)




function data_out = snicar8d_v2_dustint(BND_TYP_IN, DIRECT_IN, APRX_TYP_IN, ...
                             DELTA_IN, coszen_in, R_sfc_in, dz_in, ...
                             rho_snw_in, rds_snw_in, nbr_aer_in, ...
                             mss_cnc_sot1_in, mss_cnc_sot2_in, ...
                             mss_cnc_dst1_in, mss_cnc_dst2_in, ...
                             mss_cnc_dst3_in, mss_cnc_dst4_in, ...
                             mss_cnc_ash1_in, fl_sot1_in, fl_sot2_in, ...
                             fl_dst1_in, fl_dst2_in, fl_dst3_in, ...
                             fl_dst4_in, fl_ash1_in, ...
                             SNO_SHP_in,BC_SNO_MIX_in,SNO_fs_in,SNO_AR_in,... %che
                             mss_cnc_sot3_in,DST_SNO_MIX_in, mss_cnc_dstint_in) % che

    
if (1==0)
    % DEFINE ALL INPUT HERE (not called from external function);
    
    clear;

    % SNOW LAYER THICKNESSES [m]:
    dz = [0.02 9.98];
 
    nbr_lyr = length(dz);  % number of snow layers
  
  
    % SNOW DENSITY FOR EACH LAYER (units: kg/m3)
    rho_snw(1:nbr_lyr) = 150;  
  

    % SNOW GRAIN SIZE FOR EACH LAYER (units: microns):
    rds_snw(1:nbr_lyr) = 100;
  

    % ALTERNATIVE SNOW GRAIN SIZE: 9999
    % SET OPTICAL PROPERTIES FROM USER-DEFINED FILE:
    %rds_snw(1:nbr_lyr)=9999;
    %fl_in_usr = '/data/flanner/mie/foo/icebc_1000_0100ppb.nc';

  
    % NUMBER OF PARTICLE SPECIES IN SNOW (ICE EXCLUDED)
    nbr_aer = 7;
  
    
    % PARTICLE MASS MIXING RATIOS (units: ng g-1)
    mss_cnc_sot3(1:nbr_lyr)  = 0.0;    % che: add BC species internally mixed with snow
    mss_cnc_sot1(1:nbr_lyr)  = 0.0;    % uncoated black carbon
    mss_cnc_sot2(1:nbr_lyr)  = 0.0;    % coated black carbon
    mss_cnc_dst1(1:nbr_lyr)  = 0.0;    % dust species 1
    mss_cnc_dst2(1:nbr_lyr)  = 0.0;    % dust species 2
    mss_cnc_dst3(1:nbr_lyr)  = 0.0;    % dust species 3
    mss_cnc_dst4(1:nbr_lyr)  = 0.0;    % dust species 4
    mss_cnc_ash1(1:nbr_lyr)  = 0.0;    % volcanic ash species 1
    mss_cnc_dstint(1:nbr_lyr)  = 0.0;  % che: add dust species internally mixed with snow
  
    % FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
    % (ideally, these files should exist in all 'band' directories)
    fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
    fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
    fl_dst1  = 'aer_dst_bln_20060904_01.nc';
    fl_dst2  = 'aer_dst_bln_20060904_02.nc';
    fl_dst3  = 'aer_dst_bln_20060904_03.nc';
    fl_dst4  = 'aer_dst_bln_20060904_04.nc';
    fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';
  
    % COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM RT
    coszen=0.65;

  
    % REFLECTANCE OF SURFACE UNDERLYING SNOW:
    % (value applied to all wavelengths.  user can also specify
    % spectrally-dependent ground albedo below)
    R_sfc_all_wvl = 0;
  
  
    % RADIATIVE TRANSFER CONFIGURATION:
    BND_TYP  = 1;        % 1= 470 spectral bands, 2= 5 bands, 3= 3 bands
    DIRECT   = 1;        % 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP = 3;        % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA    = 1;        % 1= Apply Delta approximation, 0= No delta
  
    
    %%%% che %%%%%
    % Options added by Cenlin He for nonspherical snow & BC-snow internal
    % mixing based on the parameterizations from He et al. 2017 (J. Climate,doi:10.1175/JCLI-D-17-0300.1)
    
    SNO_SHP(1:nbr_lyr)  = 1;    % Snow grain shape option
                                % 1=sphere (SNICAR default shape)
                                % 2=spheroid; 3=hexagonal plate; 4=koch snowflake;

    SNO_fs(1:nbr_lyr)   = 0;    % Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                % 0=use recommended default value (He et al. 2017 JC);
                                % others(0<fs<1)= use user-specified value
                                % only activated when SNO_SHP > 1 (i.e. nonspherical)
                                
    SNO_AR(1:nbr_lyr)   = 0;    % Aspect ratio: ratio of grain width to length
                                % 0=use recommended default value (He et al. 2017 JC);
                                % others(0.1<fs<20)= use user-specified value
                                % only activated when SNO_SHP > 1 (i.e. nonspherical)
                                
    BC_SNO_MIX(1:nbr_lyr) = 1;  % BC-snow mixing option
                                % 1= only BC-snow external mixing (SNICAR default)
                                % 2= BC-snow internal + external mixing
                                % Currently, can only deal with BC-snow internal
                                % mixing; internal mixing between snow and
                                % other aerosols will be added later
    DST_SNO_MIX(1:nbr_lyr) = 1;  % dust-snow mixing option
                                % 1= only dust-snow external mixing (SNICAR default)
                                % 2= dust-snow internal + external mixing

   %%%% che %%%%%

else
    % routine has been called from external function, with needed
    % input paramaters, see snicar_dirver8d.m
    
    % RADIATIVE TRANSFER CONFIGURATION:
    BND_TYP       = BND_TYP_IN;   % 1= 470 spectral bands, 2= 5 bands, 3= 3 bands
    DIRECT        = DIRECT_IN;    % 1= Direct-beam incident flux, 0= Diffuse incident flux
    APRX_TYP      = APRX_TYP_IN;  % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
    DELTA         = DELTA_IN;     % 1= Apply Delta approximation, 0= No delta
  
    
    %%%% che %%%%% 
    % Options added by Cenlin He for nonspherical snow & BC-snow internal
    % mixing based on the parameterizations from He et al. 2017 (J. Climate,doi:10.1175/JCLI-D-17-0300.1)
    
    SNO_SHP       = SNO_SHP_in;   % 1=sphere (SNICAR default shape)
                                  % 2=spheroid; 3=hexagonal plate; 4=koch snowflake;
                                  
    SNO_fs        = SNO_fs_in;    % Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                  % 0=use recommended default value (He et al. 2017 JC);
                                  % others(0<fs<1)= use user-specified value
                                  % only activated when SNO_SHP > 1 (i.e. nonspherical)
                                  
    SNO_AR        = SNO_AR_in;    % Aspect ratio: ratio of grain width to length
                                  % 0=use recommended default value (He et al. 2017 JC);
                                  % others(0.1<fs<20)= use user-specified value
                                  % only activated when SNO_SHP > 1 (i.e. nonspherical)
                                  
    BC_SNO_MIX    = BC_SNO_MIX_in;% 1= only BC-snow external mixing (SNICAR default)
                                  % 2= BC-snow internal + external mixing
                                  % Currently, can only deal with BC-snow internal
                                  % mixing; internal mixing between snow and
                                  % other aerosols will be added later  
    DST_SNO_MIX    = DST_SNO_MIX_in;% 1= only BC-snow external mixing (SNICAR default)
                                  % 2= BC-snow internal + external mixing
    %%%% che %%%%%
    
    coszen        = coszen_in;
    R_sfc_all_wvl = R_sfc_in;

    dz            = dz_in;
    nbr_lyr       = length(dz);
    rho_snw       = rho_snw_in;
    rds_snw       = rds_snw_in;
    nbr_aer       = nbr_aer_in;

    mss_cnc_sot1  = mss_cnc_sot1_in;
    mss_cnc_sot2  = mss_cnc_sot2_in;
    mss_cnc_dst1  = mss_cnc_dst1_in;
    mss_cnc_dst2  = mss_cnc_dst2_in;
    mss_cnc_dst3  = mss_cnc_dst3_in;
    mss_cnc_dst4  = mss_cnc_dst4_in;
    mss_cnc_ash1  = mss_cnc_ash1_in;
    mss_cnc_sot3  = mss_cnc_sot3_in; %che: add BC species internally mixed with snow
    mss_cnc_dstint  = mss_cnc_dstint_in; %che: add dust species internally mixed with snow
    
    fl_sot1       = fl_sot1_in;
    fl_sot2       = fl_sot2_in;
    fl_dst1       = fl_dst1_in;
    fl_dst2       = fl_dst2_in;
    fl_dst3       = fl_dst3_in;
    fl_dst4       = fl_dst4_in;
    fl_ash1       = fl_ash1_in;

end;


% Set snow Mie directory based on band number (and direct or
% diffuse flux):
if (BND_TYP==1)
    wrkdir='./';
elseif(BND_TYP==2)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    if (DIRECT==1)
        wrkdir='/data/flanner/mie_clm/mie_clm_drc/';
    elseif (DIRECT==0)
        wrkdir='/data/flanner/mie_clm/mie_clm_dfs/';
    end;
elseif(BND_TYP==3)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    if (DIRECT==1)
        wrkdir='/data/flanner/mie_clm/mie_clm_drc_3bnd/';
    elseif (DIRECT==0)
        wrkdir='/data/flanner/mie_clm/mie_clm_dfs_3bnd/';
    end;
end;

% Set wavelength grid (um) and wavelength number:
wvl     = ncread(strcat(wrkdir,'ice_optical_properties/','ice_wrn_0100.nc'),'wvl').*10^6;
nbr_wvl = length(wvl);


% REFLECTANCE OF UNDERLYING SURFACE
% (optional: user can set spectrally-dependent surface albedo here):
R_sfc(1:nbr_wvl,1) = R_sfc_all_wvl;


% file substrings for snow/aerosol Mie parameters:
% Warren and Brandt (2008) ice optical properties:
fl_stb1 = 'ice_wrn_';
fl_stb2 = '.nc';


% Calculate mu_not
%mu_not=cos((slr_znt/360)*2*pi);
mu_not = coszen;


% Set incident solar flux spectral distribution:
if (BND_TYP==1)

    if (DIRECT == 1)
        % mid-latitude winter, clear-sky:
        load mlw_sfc_flx_frc_clr.txt;
        mlw_sfc_flx_frc_clr(find(mlw_sfc_flx_frc_clr==0))=1e-30;
        flx_slr = mlw_sfc_flx_frc_clr;
        
        % Summit Greenland, clear-sky:
        %load sas_Smm_sfc_flx_frc_clr.txt;
        %sas_Smm_sfc_flx_frc_clr(find(sas_Smm_sfc_flx_frc_clr==0))=1e-30;
        %flx_slr = sas_Smm_sfc_flx_frc_clr;

        Fs(1:nbr_wvl,1) = flx_slr./(mu_not*pi);  % direct-beam incident flux
        Fd(1:nbr_wvl,1) = 0;                     % diffuse incident flux
    
    elseif (DIRECT == 0)
        % mid-latitude winter, cloudy-sky
        load mlw_sfc_flx_frc_cld.txt;
        mlw_sfc_flx_frc_cld(find(mlw_sfc_flx_frc_cld==0))=1e-30;
        flx_slr = mlw_sfc_flx_frc_cld;
        
        % Summit Greenland, cloudy:
        %load sas_Smm_sfc_flx_frc_cld.txt;
        %sas_Smm_sfc_flx_frc_cld(find(sas_Smm_sfc_flx_frc_cld==0))=1e-30;
        %flx_slr = sas_Smm_sfc_flx_frc_cld;

        Fd(1:nbr_wvl,1) = flx_slr; % diffuse incident flux
        Fs(1:nbr_wvl,1) = 0;       % direct-beam incident flux
  end;

  
  
elseif (BND_TYP == 2)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    %
    % SPECTRAL BANDS:
    %  1. 0.3-0.7 um
    %  2. 0.7-1.0 um
    %  3. 1.0-1.2 um
    %  4. 1.2-1.5 um
    %  5. 1.5-5.0 um

    % mid-latitude winter:
    if (DIRECT == 1)
        flx_slr(1,1) = 0.5028;
        flx_slr(2,1) = 0.2454;
        flx_slr(3,1) = 0.0900;
        flx_slr(4,1) = 0.0601;
        flx_slr(5,1) = 0.1017;
   
        Fs = flx_slr./(mu_not*pi);
        Fd(1:nbr_wvl,1) = 0;
  
    elseif (DIRECT == 0)
        flx_slr(1,1) = 0.5767;
        flx_slr(2,1) = 0.2480;
        flx_slr(3,1) = 0.0853;
        flx_slr(4,1) = 0.0462;
        flx_slr(5,1) = 0.0438;
        
        Fs(1:nbr_wvl,1) = 0;
        Fd = flx_slr;
    end;
    
    
elseif (BND_TYP == 3)
    % NOT FUNCTIONAL IN THIS DISTRIBUTION
    %
    % SPECTRAL BANDS:
    %  1. 0.3-0.7 um
    %  2. 0.7-1.19 um
    %  3. 1.0-5.0 um
      
    if (DIRECT == 1)
        % MID-LATITUDE WINTER (SWNB2) (0.3-0.7, 0.7-1.19, 1.19-5.0)
        flx_slr(1,1) = 0.5028;
        flx_slr(2,1) = 0.3313;
        flx_slr(3,1) = 0.1659;
   
        Fs = flx_slr./(mu_not*pi);
        Fd(1:nbr_wvl,1) = 0;
  
    elseif (DIRECT == 0)
        flx_slr(1,1) = 0.5767;
        flx_slr(2,1) = 0.3297;
        flx_slr(3,1) = 0.0936;
        
        Fs(1:nbr_wvl,1) = 0;
        Fd = flx_slr;
    end;
    
end;


% Read Mie ice parameters from external NetCDF files (layer-specific):
ext_cff_mss_snw = zeros(nbr_wvl, nbr_lyr);
omega_snw = zeros(nbr_wvl, nbr_lyr);
g_snw = zeros(nbr_wvl, nbr_lyr);
R_1_omega = zeros(nbr_wvl,nbr_lyr); % BC-induced enhancement in snow single-scattering coalbedo (1-omega)


%%% che: Added by Cenlin He %%%
% g_snw asymmetry factor parameterization coefficients (6 bands) from
% Table 3 & Eqs. 6-7 in He et al. 2017 J. Climate
% assume same values for 4-5 um band, which leads to very small biases (<3%)
g_wvl = [0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00]; % wavelength (um) division point
g_wvl_center = g_wvl(2:8)/2 + g_wvl(1:7)/2 ; % center point for wavelength band
g_b0 = [9.76029E-01,9.67798E-01,1.00111E+00,1.00224E+00,9.64295E-01,9.97475E-01,9.97475E-01];
g_b1 = [5.21042E-01,4.96181E-01,1.83711E-01,1.37082E-01,5.50598E-02,8.48743E-02,8.48743E-02];
g_b2 = [-2.66792E-04,1.14088E-03,2.37011E-04,-2.35905E-04,8.40449E-04,-4.71484E-04,-4.71484E-04];
% Tables 1 & 2 and Eqs. 3.1-3.4 in Fu, 2007 
g_F07_c2 = [1.349959e-1,1.115697e-1,9.853958e-2,5.557793e-2,-1.233493e-1,0.0,0.0];
g_F07_c1 = [-3.987320e-1,-3.723287e-1,-3.924784e-1,-3.259404e-1,4.429054e-2,-1.726586e-1,-1.726586e-1];
g_F07_c0 = [7.938904e-1,8.030084e-1,8.513932e-1,8.692241e-1,7.085850e-1,6.412701e-1,6.412701e-1];
g_F07_p2 = [3.165543e-3,2.014810e-3,1.780838e-3,6.987734e-4,-1.882932e-2,-2.277872e-2,-2.277872e-2];
g_F07_p1 = [1.140557e-1,1.143152e-1,1.143814e-1,1.071238e-1,1.353873e-1,1.914431e-1,1.914431e-1];
g_F07_p0 = [5.292852e-1,5.425909e-1,5.601598e-1,6.023407e-1,6.473899e-1,4.634944e-1,4.634944e-1];
			    
% Parameterization coefficients for BC-induced enhancement in snow single-scattering coalbedo
% Eq. 8b & Table 4 in He et al., 2017 J. Climate
% lambda>1.2um, no BC enhancement, enhancement ratio = 1.0
R_wvl = [0.25,0.30,0.33,0.36,0.40,0.44,0.48,0.52,0.57,0.64,0.69,0.75,0.78,...
    0.87,1.00,1.20,5.00]; % wavelength (um) division point
R_d0 = [4.70305E+00,4.68619E+00,4.67369E+00,4.65040E+00,2.40364E+00,...
    7.95408E-01,2.92745E-01,8.63396E-02,2.76299E-02,1.40864E-02,...
    8.65705E-03,6.12971E-03,4.45697E-03,3.06648E-02,7.96544E-01,1.0];
R_d1 = [9.73317E-01,9.79650E-01,9.84579E-01,9.93537E-01,9.95955E-01,...
    9.95218E-01,9.74284E-01,9.81193E-01,9.81239E-01,9.55515E-01,...
    9.10491E-01,8.74196E-01,8.27238E-01,4.82870E-01,4.36649E-02,0.0];
R_d2 = [2.04820E-01,2.07410E-01,2.09390E-01,2.13030E-01,4.18570E-01,...
    1.29682E+00,3.75514E+00,1.27372E+01,3.93293E+01,8.78918E+01,...
    1.86969E+02,3.45600E+02,7.08637E+02,1.41067E+03,2.57288E+02,1.0];
nwvl_R = length(R_d0);
R_wvl_center = R_wvl(2:(nwvl_R+1))/2 + R_wvl(1:nwvl_R)/2 ; % center point for wavelength band

% Parameterization coefficients for dust-induced enhancement in snow single-scattering coalbedo
% Eq. 1 & Table 1 in He et al., 2019 JAMES
% lambda>1.2um, no dust enhancement, enhancement ratio = 1.0
E_wvl = [0.2000,0.2632,0.3448,0.4415,0.6250,0.7782,1.2422];
E_a1  = [-2.1307E+01,-1.5815E+01,-9.2880E+00,1.1115E+00,1.0307E+00,1.0185E+00] ;
E_a2  = [ 1.1746E+02, 9.3241E+01, 4.0605E+01,3.7389E-01,1.4800E-02,2.8921E-04] ;
E_a3  = [ 9.9701E-01, 9.9781E-01, 9.9848E-01,1.0035E+00,1.0024E+00,1.0356E+00] ;
nwvl_E = length(E_a1);
E_wvl_center = E_wvl(2:(nwvl_E+1))/2 + R_wvl(1:nwvl_E)/2 ; % center point for wavelength band
%%%


for n=1:nbr_lyr
    if (rds_snw(n)<10)
        s1    = int2str(0);
        s2    = int2str(rds_snw(n));
        fl_in = strcat(wrkdir,'ice_optical_properties/',fl_stb1,s1,s1,s1,s2,fl_stb2);
    elseif (rds_snw(n)<100)
        s1    = int2str(0);
        s2    = int2str(rds_snw(n));
        fl_in = strcat(wrkdir,'ice_optical_properties/',fl_stb1,s1,s1,s2,fl_stb2);
    elseif ((rds_snw(n)>=100) & (rds_snw(n)<1000))
        s1    = int2str(0);
        s2    = int2str(rds_snw(n));
        fl_in = strcat(wrkdir,'ice_optical_properties/',fl_stb1,s1,s2,fl_stb2);
    elseif ((rds_snw(n)>=1000) & (rds_snw(n)<9999))
        s1    = int2str(0);
        s2    = int2str(rds_snw(n));
        fl_in = strcat(wrkdir,'ice_optical_properties/',fl_stb1,s2,fl_stb2);
    elseif (rds_snw(n)==9999)
        fl_in = fl_in_usr;
    end;
    
    % mass extinction coefficient & single-scatter albedo not affected by snow shape
    ext_cff_mss_snw(:,n) = ncread(fl_in,'ext_cff_mss');
    omega_snw(:,n)       = ncread(fl_in,'ss_alb'); 
    
    %%% che: added by Cenlin He for asymmetry factor affected by snow shapes %%%%
    if (SNO_SHP(n) == 1)  % SNICAR default sphere from MIE look-up table
        g_snw(:,n)       = ncread(fl_in,'asm_prm'); % asymmetry paramater
        
    elseif (SNO_SHP(n) == 2) % 2=spheroid, He et al. 2017 parameterization
        diam_snw = 2.0 .* rds_snw(n); % effective snow grain diameter        
        if (SNO_fs(n) == 0)
            fs_sphd = 0.929; % default shape factor for spheroid; He et al. 2017 Table 1
        else
            fs_sphd = SNO_fs(n);
        end
        fs_hex = 0.788; % shape factor for hexagonal plate (reference);
        if (SNO_AR(n) == 0)
            AR_tmp = 0.5; % default aspect ratio for spheroid; He et al. 2017 Table 1
        else
            AR_tmp = SNO_AR(n);
        end
        g_snw_Cg_tmp = g_b0 .* (fs_sphd/fs_hex).^g_b1 .* diam_snw.^g_b2; % Eq.7, He et al. 2017
        gg_snw_F07_tmp = g_F07_c0 + g_F07_c1 .* AR_tmp + g_F07_c2 .* AR_tmp^2;% Eqs.3.1 in Fu, 2007
        % 6 wavelength bands for g_snw to be interpolated into 470-bands of SNICAR
        % shape-preserving piecewise interpolation into 470-bands 
        g_Cg_intp = pchip(g_wvl_center,g_snw_Cg_tmp,wvl) ;
        gg_F07_intp = pchip(g_wvl_center,gg_snw_F07_tmp,wvl) ;
        g_snw_F07 = gg_F07_intp + (1.0 - gg_F07_intp) ./ omega_snw(:,n) ./ 2; % Eq.2.2 in Fu, 2007
        g_snw(:,n) = g_snw_F07 .* g_Cg_intp; % Eq.6, He et al. 2017
        g_snw(371:470,n) = g_snw(370); % assume same values for 4-5 um band, with very small biases (<3%)
        
    elseif (SNO_SHP(n) == 3) % 3=hexagonal plate, He et al. 2017 parameterization
        diam_snw = 2.0 .* rds_snw(n); % effective snow grain diameter        
        if (SNO_fs(n) == 0)
            fs_hex0 = 0.788; % default shape factor for hexagonal plates; He et al. 2017 Table 1
        else
            fs_hex0 = SNO_fs(n);
        end
        fs_hex = 0.788; % shape factor for hexagonal plate (reference);
        if (SNO_AR(n) == 0)
            AR_tmp = 2.5; % default aspect ratio for hexagonal plate; He et al. 2017 Table 1
        else
            AR_tmp = SNO_AR(n);
        end
        g_snw_Cg_tmp = g_b0 .* (fs_hex0/fs_hex).^g_b1 .* diam_snw.^g_b2; % Eq.7, He et al. 2017
        gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 .* log(AR_tmp) + g_F07_p2 .* (log(AR_tmp))^2;   % Eqs.3.3 in Fu, 2007 
        % 6 wavelength bands for g_snw to be interpolated into 470-bands of SNICAR
        % shape-preserving piecewise interpolation into 470-bands 
        g_Cg_intp = pchip(g_wvl_center,g_snw_Cg_tmp,wvl) ;
        gg_F07_intp = pchip(g_wvl_center,gg_snw_F07_tmp,wvl) ;
        g_snw_F07 = gg_F07_intp + (1.0 - gg_F07_intp) ./ omega_snw(:,n) ./ 2; % Eq.2.2 in Fu, 2007
        g_snw(:,n) = g_snw_F07 .* g_Cg_intp; % Eq.6, He et al. 2017
        g_snw(371:470,n) = g_snw(370); % assume same values for 4-5 um band, with very small biases (<3%)
        
    elseif (SNO_SHP(n) == 4) % 4=koch snowflake, He et al. 2017 parameterization
        diam_snw = 2.0 .* rds_snw(n) ./ 0.544; % effective snow grain diameter
        if (SNO_fs(n) == 0)
            fs_koch = 0.712; % default shape factor for koch snowflake; He et al. 2017 Table 1
        else
            fs_koch = SNO_fs(n);
        end
        fs_hex = 0.788; % shape factor for hexagonal plate (reference);
        if (SNO_AR(n) == 0)
            AR_tmp = 2.5; % default aspect ratio for koch snowflake; He et al. 2017 Table 1
        else
            AR_tmp = SNO_AR(n);
        end
        g_snw_Cg_tmp = g_b0 .* (fs_koch/fs_hex).^g_b1 .* diam_snw.^g_b2; % Eq.7, He et al. 2017
        gg_snw_F07_tmp = g_F07_p0 + g_F07_p1 .* log(AR_tmp) + g_F07_p2 .* (log(AR_tmp))^2;   % Eqs.3.3 in Fu, 2007 
        % 6 wavelength bands for g_snw to be interpolated into 470-bands of SNICAR
        % shape-preserving piecewise interpolation into 470-bands 
        g_Cg_intp = pchip(g_wvl_center,g_snw_Cg_tmp,wvl) ;
        gg_F07_intp = pchip(g_wvl_center,gg_snw_F07_tmp,wvl) ;
        g_snw_F07 = gg_F07_intp + (1.0 - gg_F07_intp) ./ omega_snw(:,n) ./ 2; % Eq.2.2 in Fu, 2007
        g_snw(:,n) = g_snw_F07 .* g_Cg_intp; % Eq.6, He et al. 2017
        g_snw(371:470,n) = g_snw(370); % assume same values for 4-5 um band, with very small biases (<3%)
        
    end
    
    g_snw(g_snw > 0.99) = 0.99; % avoid unreasonable values (so far only occur in large-size spheroid cases)
    
    %%% che: added by Cenlin He for BC-snow internal mixing
    if (BC_SNO_MIX(n) == 2 && mss_cnc_sot3(n) > 0)   % BC-snow internal mixing
        % R_1_omega: BC-induced enhancement in snow single-scattering coalbedo
        % R_1_omega from Eq.8b in He et al.(2017,JC) is based on BC Re=0.1um &
        % MAC=6.81 m2/g (@550 nm) & BC density=1.7g/cm3.
        % To be consistent with SNICAR default (BC MAC=7.5 m2/g @550nm), we
        % made adjustments on BC size & density as follows to get MAC=7.5m2/g.
        % (1) We use BC Re=0.045um [geometric mean diameter=0.06um (Dentener et al.2006, 
        % Yu and Luo,2009) & geometric std=1.5 (Flanner et al.2007;Aoki et al., 2011)]
        % (2) We tune BC density from 1.7 to 1.49 g/cm3 (Aoki et al., 2011) to match BC MAC=7.5 m2/g @550 nm.
        mss_cnc_sot3(n) = mss_cnc_sot3(n) * 1.7/1.49; % BC density adjustment
        R_1_omega_tmp = R_d0 .* (mss_cnc_sot3(n) + R_d2).^R_d1 ; % Eq. 8b in He et al.2017,JC
        % Adust R_1_omega_tmp due to BC Re from 0.1 to 0.045um based on 
        % Eq. 1 & Table S1 in He et al.2018 (GRL)
        R_1_omega_010_vis = R_1_omega_tmp(1:10); % enhancement for Re=0.1um, visible
        R_1_omega_010_nir = R_1_omega_tmp(11:15); % enhancement for Re=0.1um, NIR
        R_1_omega_005_vis = (R_1_omega_010_vis ./ (0.1/0.05)^-0.1866).^((0.1/0.05)^0.1918) ; % visible
        R_1_omega_005_nir = (R_1_omega_010_nir ./ (0.1/0.05)^-0.0046).^((0.1/0.05)^0.5177) ; % NIR
        R_1_omega_0045_vis = (0.045/0.05)^-0.1866 .* (R_1_omega_005_vis.^((0.045/0.05)^-0.1918)); % visible
        R_1_omega_0045_nir = (0.045/0.05)^-0.0046 .* (R_1_omega_005_nir.^((0.045/0.05)^-0.5177)); % NIR
        R_1_omega_tmp(1:10) = R_1_omega_0045_vis; % visible
        R_1_omega_tmp(11:15) = R_1_omega_0045_nir; % NIR
        % shape-preserving piecewise interpolation into 470-bands
        logR_1_omega_tmp = log10(R_1_omega_tmp); % interpolate in a logscale field
        logR_intp = pchip(R_wvl_center,logR_1_omega_tmp,wvl) ;
        R_intp = 10.0 .^ logR_intp;
        R_intp(R_intp < 1.0) = 1.0; % BC does not reduce snow absorption       
        R_intp(R_intp > 1e5) = 1e5;       
        R_intp(91:470) = 1.0; % lambda>1.2um, no BC enhancement
        R_1_omega(:,n) = R_intp;
        % new 1-omega (single-scattering coalbedo) for entire BC-snow internal mixture
        tmp_1_omega_snw = (1.0 - omega_snw(:,n)) .* R_1_omega(:,n); 
        % Recommended: smoothing the enhanced 1-omega for <1.2um (small
        % wiggling pattern would occur sometime). Smoothing has negligible 
        % effects on broadband (vis & NIR) omega & albedo
     %   tmp_1_omega_sm = smooth(log10(tmp_1_omega_snw(1:90)),5) ; % 5-point moving average smooth
     %   tmp_1_omega_snw(1:90) = 10.0 .^ tmp_1_omega_sm ; 
        tmp_data = tmp_1_omega_snw;
        tmp_data(91:470) = 0.0; % do not treat wavelength>1.2um
        tmp_1_omega_snw(tmp_data>tmp_1_omega_snw(91)) = tmp_1_omega_snw(91); % prevenet snow ssa to go too low and keep continuity
        % new omega for entire BC-snow internal mixture
        omega_snw(:,n) = 1.0 - tmp_1_omega_snw;
    end
    %%%

    %%% che: added by Cenlin He for dust-snow internal mixing
    if (DST_SNO_MIX(n) == 2 && mss_cnc_dstint(n) > 0)   % dust-snow internal mixing
        % E_1_omega: dust-induced enhancement in snow single-scattering coalbedo
        % E_1_omega from Eq.1 in He et al.(2019,JAMES), Cdust: ppm unit
        E_1_omega_tmp = E_a1 + E_a2 .* (mss_cnc_dstint(n)/1000).^E_a3 ; 
        % shape-preserving piecewise interpolation into 470-bands
        logE_1_omega_tmp = log10(E_1_omega_tmp); % interpolate in a logscale field
        logE_intp = pchip(E_wvl_center,logE_1_omega_tmp,wvl) ;
        E_intp = 10.0 .^ logE_intp;
        E_intp(E_intp < 1.0) = 1.0; % dust does not reduce snow absorption       
        E_intp(E_intp > 1e5) = 1e5;        
        E_intp(91:470) = 1.0; % lambda>1.2um, no dust enhancement
        % new 1-omega (single-scattering coalbedo) for entire dust-snow internal mixture
        tmp_1_omega_snw = (1.0 - omega_snw(:,n)) .* E_intp;
        tmp_data = tmp_1_omega_snw;
        tmp_data(91:470) = 0.0; % do not treat wavelength>1.2um
        tmp_1_omega_snw(tmp_data>tmp_1_omega_snw(91)) = tmp_1_omega_snw(91); % prevenet snow ssa to go too low and keep continuity
        % Recommended: smoothing the enhanced 1-omega for <1.2um (small
        % wiggling pattern would occur sometime). Smoothing has negligible 
        % effects on broadband (vis & NIR) omega & albedo
      %  tmp_1_omega_sm = smooth(log10(tmp_1_omega_snw(1:90)),5) ; % 5-point moving average smooth
      %  tmp_1_omega_snw(1:90) = 10.0 .^ tmp_1_omega_sm ; 
        % new omega for entire BC-snow internal mixture
        omega_snw(:,n) = 1.0 - tmp_1_omega_snw;
    end
    %%%
 
end

% Read Mie aerosol parameters (layer-independent)
fl_in1 = strcat(wrkdir,fl_sot1);
fl_in2 = strcat(wrkdir,fl_sot2);
fl_in3 = strcat(wrkdir,fl_dst1);
fl_in4 = strcat(wrkdir,fl_dst2);
fl_in5 = strcat(wrkdir,fl_dst3);
fl_in6 = strcat(wrkdir,fl_dst4);
fl_in7 = strcat(wrkdir,fl_ash1);

omega_aer(:,1)       = ncread(fl_in1,'ss_alb');
g_aer(:,1)           = ncread(fl_in1,'asm_prm');
ext_cff_mss_aer(:,1) = ncread(fl_in1,'ext_cff_mss');
  
omega_aer(:,2)       = ncread(fl_in2,'ss_alb');
g_aer(:,2)           = ncread(fl_in2,'asm_prm');
ext_cff_mss_aer(:,2) = ncread(fl_in2,'ext_cff_mss_cor');
%ext_cff_mss_aer(:,2)= ncread(fl_in2,'ext_cff_mss');

omega_aer(:,3)       = ncread(fl_in3,'ss_alb');
g_aer(:,3)           = ncread(fl_in3,'asm_prm');
ext_cff_mss_aer(:,3) = ncread(fl_in3,'ext_cff_mss');

omega_aer(:,4)       = ncread(fl_in4,'ss_alb');
g_aer(:,4)           = ncread(fl_in4,'asm_prm');
ext_cff_mss_aer(:,4) = ncread(fl_in4,'ext_cff_mss');

omega_aer(:,5)       = ncread(fl_in5,'ss_alb');
g_aer(:,5)           = ncread(fl_in5,'asm_prm');
ext_cff_mss_aer(:,5) = ncread(fl_in5,'ext_cff_mss');

omega_aer(:,6)       = ncread(fl_in6,'ss_alb');
g_aer(:,6)           = ncread(fl_in6,'asm_prm');
ext_cff_mss_aer(:,6) = ncread(fl_in6,'ext_cff_mss');

omega_aer(:,7)       = ncread(fl_in7,'ss_alb');
g_aer(:,7)           = ncread(fl_in7,'asm_prm');
ext_cff_mss_aer(:,7) = ncread(fl_in7,'ext_cff_mss');


% Set aerosol concentration matrix:
mss_cnc_aer(1:nbr_lyr,1) = mss_cnc_sot1;
mss_cnc_aer(1:nbr_lyr,2) = mss_cnc_sot2;
mss_cnc_aer(1:nbr_lyr,3) = mss_cnc_dst1;
mss_cnc_aer(1:nbr_lyr,4) = mss_cnc_dst2;
mss_cnc_aer(1:nbr_lyr,5) = mss_cnc_dst3;
mss_cnc_aer(1:nbr_lyr,6) = mss_cnc_dst4;
mss_cnc_aer(1:nbr_lyr,7) = mss_cnc_ash1;
mss_cnc_aer(1:nbr_lyr,7) = mss_cnc_ash1;


% convert to units of kg/kg:
mss_cnc_aer = mss_cnc_aer.*10^-9;


% BEGIN RT SOLVER:

% Calculate effective tau, omega, g for the (snow+impurity) system
for n=1:nbr_lyr
    L_snw(n)     = rho_snw(n)*dz(n);                 % Snow column mass (kg/m^2) (array)
    tau_snw(:,n) = L_snw(n).*ext_cff_mss_snw(:,n);
    
    for j=1:nbr_aer
        L_aer(n,j)     = L_snw(n)*mss_cnc_aer(n,j);
        tau_aer(:,n,j) = L_aer(n,j).*ext_cff_mss_aer(:,j);
    end
  
    tau_sum(1:nbr_wvl,1)   = 0.0;
    omega_sum(1:nbr_wvl,1) = 0.0;
    g_sum(1:nbr_wvl,1)     = 0.0;
  
    for j=1:nbr_aer
        tau_sum   = tau_sum + tau_aer(:,n,j);
        omega_sum = omega_sum + (tau_aer(:,n,j).*omega_aer(:,j));
        g_sum     = g_sum + (tau_aer(:,n,j).*omega_aer(:,j).*g_aer(:,j));
    end
    
    tau(:,n)   = tau_sum + tau_snw(:,n);
    omega(:,n) = (1./tau(:,n)).*(omega_sum+ (omega_snw(:,n).*tau_snw(:,n)));
    g(:,n)     = (1./(tau(:,n).*omega(:,n))) .* (g_sum+ (g_snw(:,n).*omega_snw(:,n).*tau_snw(:,n)));
end
  
  
% Perform Delta-transformations, if called for
if (DELTA == 1)
    g_star     = g./(1+g);
    omega_star = ((1-(g.^2)).*omega) ./ (1-(omega.*(g.^2)));
    tau_star   = (1-(omega.*(g.^2))).*tau;
else
    g_star     = g;
    omega_star = omega;
    tau_star   = tau;
end;

% Calculate total column optical depth
% tau_clm(:,n) = total optical depth of layers above n, or optical
% depth from upper model boundary to top of layer n
tau_clm(1:nbr_wvl,1:1) = 0.0;
for n=2:nbr_lyr
    tau_clm(:,n) = tau_clm(:,n-1)+tau_star(:,n-1);
end


% Boundary Condition: Upward (reflected) direct beam flux at lower model boundary
% (Eq. 37)
S_sfc = R_sfc.*mu_not.*exp(-(tau_clm(:,nbr_lyr)+tau_star(:,nbr_lyr))./mu_not).*pi.*Fs;


% Apply 2-stream approximation technique (Toon et al., Table 1.)

if (APRX_TYP == 1)
    % Eddington:
    gamma1 = (7-(omega_star.*(4+(3*g_star))))./4;
    gamma2 = -(1-(omega_star.*(4-(3*g_star))))./4;
    gamma3 = (2-(3*g_star.*mu_not))./4;
    gamma4 = 1-gamma3;
    mu_one = 0.5;
  
elseif (APRX_TYP == 2)
    % Quadrature: 
    gamma1 = sqrt(3).*(2-(omega_star.*(1+g_star)))./2;
    gamma2 = omega_star.*sqrt(3).*(1-g_star)./2;
    gamma3 = (1-(sqrt(3).*g_star.*mu_not))./2;
    gamma4 = 1-gamma3;
    mu_one = 1/sqrt(3);
  
elseif (APRX_TYP == 3)
    % Hemispheric mean:
    gamma1 = 2 - (omega_star.*(1+g_star));
    gamma2 = omega_star.*(1-g_star);
    gamma3 = (1-(sqrt(3).*g_star.*mu_not))./2;
    gamma4 = 1-gamma3;
    mu_one = 0.5;
end;

% Eq. 21 and 22
lambda = sqrt(abs((gamma1.^2) - (gamma2.^2)));
GAMMA  = gamma2./(gamma1+lambda);

% Eq. 44
e1 = 1+(GAMMA.*exp(-lambda.*tau_star));
e2 = 1-(GAMMA.*exp(-lambda.*tau_star));
e3 = GAMMA + exp(-lambda.*tau_star);
e4 = GAMMA - exp(-lambda.*tau_star);


% Before calculating the C functions:
% C functions are indeterminate if [lambda^2 = 1/(mu_not^2)]
% upper bound of lambda^2 = 4 for delta-hemispheric mean, and
% lambda^2 = 3 for delta-Eddington.  So, problem can only arise
% when 0.5 < mu_not < 1.0.  
% Assuming that abs(lambda^2 - 1/(mu_not^2)) < 0.01 is dangerous:
% Let x= lambda^2 - 1/(mu_not^2),
% dx/dmu_not = 2*mu_not^-3
% dx = 0.01 = 2*(mu_not^-3)*dmu_not
% For range of possible mu_not:
% 0.04 > dmu_not > 0.005
% So, changing mu_not by 0.04 will be sufficient:

% Note: not implemented here because of vectorized code, 
% but an example is shown.  MUST be implemented in CLM

% Calculate C function for each layer
% (evaluated at the top and bottom of each layer n)
% (Eq. 23 and 24)

for n=1:nbr_lyr

    if (sum(Fs) > 0.0)
      
        % Stability check, see explanation above:
        %for p=1:nbr_wvl
        %temp1 = ( lambda(p,n)^2) - 1/mu_not^2);
        %if abs(temp) < 0.01
        %  mu_not_sv = mu_not;
        %  if (temp1 > 0)
        %    mu_not = mu_not + 0.04;
        %  else
        %    mu_not = mu_not - 0.04
        %  end;
        %end;
        %end
        %if (n==1)
        %    mu_not=mu_not+0.00;
        %else
        %    mu_not = mu_not-0.00;
        %end;
	
        C_pls_btm(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not).*...
                          (((gamma1(:,n)-(1/mu_not)).*gamma3(:,n))+ ...
                           (gamma4(:,n).*gamma2(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
  
            
        C_mns_btm(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not).*...
                          (((gamma1(:,n)+(1/mu_not)).*gamma4(:,n))+ ...
                           (gamma2(:,n).*gamma3(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));


        C_pls_top(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-tau_clm(:,n)./mu_not).*...
                          (((gamma1(:,n)-(1/mu_not)).*gamma3(:,n))+(gamma4(:,n).*gamma2(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
        
            
        C_mns_top(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-tau_clm(:,n)./mu_not).*...
                          (((gamma1(:,n)+(1/mu_not)).*gamma4(:,n))+(gamma2(:,n).*gamma3(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
        
    else
        % no direct-beam flux:
        C_pls_btm(1:nbr_wvl,n) = 0.0;
        C_mns_btm(1:nbr_wvl,n) = 0.0;
        C_pls_top(1:nbr_wvl,n) = 0.0;
        C_mns_top(1:nbr_wvl,n) = 0.0;
    end;
      
end
  

% Eq. 41-43
for i=1:2*nbr_lyr
    % Boundary values for i=1 and i=2nbr_lyr, specifics for i=odd and i=even    
    if (i==1)
        A(1:nbr_wvl,1) = 0.0;
        B(1:nbr_wvl,1) = e1(:,1);
        D(1:nbr_wvl,1) = -e2(:,1);
        E(1:nbr_wvl,1) = Fd(:)-C_mns_top(:,1);
					  
    elseif(i==2*nbr_lyr)
        A(:,i) = e1(:,nbr_lyr)-(R_sfc.*e3(:,nbr_lyr));
        B(:,i) = e2(:,nbr_lyr)-(R_sfc.*e4(:,nbr_lyr));
        D(:,i) = 0.0;
        E(:,i) = S_sfc(:) - C_pls_btm(:,nbr_lyr) + (R_sfc.*C_mns_btm(:,nbr_lyr));
    
    elseif(mod(i,2)==1)  % If odd and i>=3 (n=1 for i=3. this is confusing)
        n      = floor(i/2);
        A(:,i) = (e2(:,n).*e3(:,n))-(e4(:,n).*e1(:,n));
        B(:,i) = (e1(:,n).*e1(:,n+1))-(e3(:,n).*e3(:,n+1));
        D(:,i) = (e3(:,n).*e4(:,n+1))-(e1(:,n).*e2(:,n+1));
        E(:,i) = (e3(:,n).*(C_pls_top(:,n+1)-C_pls_btm(:,n))) + ...
                 (e1(:,n).*(C_mns_btm(:,n)-C_mns_top(:,n+1)));
    
    elseif(mod(i,2)==0)  % If even and i<=2nbr_lyr
        n      = (i/2);
        A(:,i) = (e2(:,n+1).*e1(:,n))-(e3(:,n).*e4(:,n+1));
        B(:,i) = (e2(:,n).*e2(:,n+1))-(e4(:,n).*e4(:,n+1));
        D(:,i) = (e1(:,n+1).*e4(:,n+1))-(e2(:,n+1).*e3(:,n+1));
        E(:,i) = (e2(:,n+1).*(C_pls_top(:,n+1)-C_pls_btm(:,n))) + ...
                 (e4(:,n+1).*(C_mns_top(:,n+1)-C_mns_btm(:,n))); 
    end;
end


% Eq. 45
AS(1:nbr_wvl,2*nbr_lyr) = A(:,2*nbr_lyr)./B(:,2*nbr_lyr);
DS(1:nbr_wvl,2*nbr_lyr) = E(:,2*nbr_lyr)./B(:,2*nbr_lyr);
  
% Eq. 46
for i=(2*nbr_lyr-1):-1:1
    X(1:nbr_wvl,i) = 1./(B(:,i)-(D(:,i).*AS(:,i+1)));
    AS(:,i)        = A(:,i).*X(:,i);
    DS(:,i)        = (E(:,i)-(D(:,i).*DS(:,i+1))).*X(:,i);
end

% Eq. 47
Y(1:nbr_wvl,1) = DS(:,1);
for i=2:2*nbr_lyr
    Y(:,i) = DS(:,i)-(AS(:,i).*Y(:,i-1));
end;


for n=1:nbr_lyr
    % Direct beam flux at the base of each layer (Eq. 50)
    direct(1:nbr_wvl,n) = mu_not*pi*Fs.*exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not);

    % Net flux (positive upward = F_up-F_down) at the base of each
    % layer (Eq. 48)
    F_net(1:nbr_wvl,n) = (Y(:,(2*n-1)).*(e1(:,n)-e3(:,n))) +...
        (Y(:,(2*n)).*(e2(:,n)-e4(:,n))) + ...
        C_pls_btm(:,n) - C_mns_btm(:,n) - direct(:,n);

    % Mean intensity at the base of each layer (Eq. 49):
    intensity(1:nbr_wvl,n) = (1/mu_one).*...
        ( Y(:,(2*n-1)).*(e1(:,n)+e3(:,n)) + ...
          Y(:,(2*n)).*(e2(:,n)+e4(:,n)) + C_pls_btm(:,n) + C_mns_btm(:,n)) +...
        (direct(:,n)./mu_not);
  
    intensity(1:nbr_wvl,n) = intensity(1:nbr_wvl,n)./(4*pi);
end
  

% Upward flux at upper model boundary (Eq. 31):
F_top_pls = (Y(:,1).*(exp(-lambda(:,1).*tau_star(:,1))+GAMMA(:,1))) + ...
    (Y(:,2).*(exp(-lambda(:,1).*tau_star(:,1))-GAMMA(:,1))) + ...
    C_pls_top(:,1);


for n=1:nbr_lyr
    % Upward flux at the bottom of each layer interface (Eq. 31)
    F_up(1:nbr_wvl,n) = ...
        Y(:,2*n-1).*(exp(0) + GAMMA(:,n).*exp(-lambda(:,n).*tau_star(:,n))) +...
        Y(:,2*n)  .*(exp(0) - GAMMA(:,n).*exp(-lambda(:,n).*tau_star(:,n))) +...
        C_pls_btm(:,n);
  
    % Downward flux at the bottom of each layer interface (Eq. 32,
    % plus direct-beam component):
    F_down(1:nbr_wvl,n) = ...
        Y(:,2*n-1).*(GAMMA(:,n).*exp(0) + exp(-lambda(:,n).*tau_star(:,n))) + ...
        Y(:,2*n)  .*(GAMMA(:,n).*exp(0) - exp(-lambda(:,n).*tau_star(:,n))) + ...
        C_mns_btm(:,n) + direct(:,n);

    % Derived net (upward-downward) flux
    % (should equal F_net, Eq. 48)
    F_net2(1:nbr_wvl,n) = F_up(1:nbr_wvl,n) - F_down(1:nbr_wvl,n);
  
    % planar intensity:
    intensity2(1:nbr_wvl,n) = F_up(1:nbr_wvl,n) + F_down(1:nbr_wvl,n);
end

% surface planar intensity:
intensity2_top(1:nbr_wvl) = F_top_pls + ((mu_not*pi*Fs)+Fd);

% diagnostic:
if (BND_TYP==1)
    intensity_out(1)           = intensity2_top(26);
    intensity_out(2:nbr_lyr+1) = intensity2(26,1:nbr_lyr);
end;


% Net flux at lower model boundary = bulk transmission through entire
% media = absorbed radiation by underlying surface:
F_btm_net = -F_net(:,nbr_lyr);


% Hemispheric wavelength-dependent albedo:
if (BND_TYP < 4)
    albedo = F_top_pls./((mu_not*pi*Fs)+Fd);
end


% Net flux at upper model boundary
F_top_net(1:nbr_wvl,1) = F_top_pls - ((mu_not*pi*Fs)+Fd);

  
% Absorbed flux in each layer (negative if there is net emission (bnd_typ==4))
for n=1:nbr_lyr
    if(n==1)
        F_abs(1:nbr_wvl,1) = F_net(:,1)-F_top_net;
    else
        F_abs(:,n) = F_net(:,n) - F_net(:,n-1);
    end;
end


% Set indices for VIS and NIR
if (BND_TYP == 1)
    vis_max_idx = 40;
    nir_max_idx = length(wvl);
elseif ((BND_TYP == 2) | (BND_TYP == 3))
    vis_max_idx = 1;
    nir_max_idx = length(wvl);
end;


% Spectrally-integrated absorption in each layer:
F_abs_slr=sum(F_abs);
for n=1:nbr_lyr
    F_abs_vis(n) = sum(F_abs(1:vis_max_idx,n));
    F_abs_nir(n) = sum(F_abs(vis_max_idx+1:nir_max_idx,n));
end

% Spectrally-integrated absorption by underlying surface:
F_abs_btm = sum(F_btm_net);
F_abs_vis_btm = sum(F_btm_net(1:vis_max_idx));
F_abs_nir_btm = sum(F_btm_net(vis_max_idx+1:nir_max_idx));


% Radiative heating rate:
heat_rt = F_abs_slr./(L_snw.*2117);   % [K/s] 2117 = specific heat ice (J kg-1 K-1)	   
heat_rt = heat_rt.*3600;              %[K/hr]
				      
				      
% Energy conservation check:
% Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
energy_sum = (mu_not*pi*Fs)+Fd - (sum(F_abs,2) + F_btm_net + F_top_pls);

if (sum(abs(energy_sum)) > 1e-10)
    energy_conservation_error = sum(abs(energy_sum))
    %error(strcat('Energy conservation error of: ',num2str(sum(abs(energy_sum)))));
end;

% spectrally-integrated terms (remove semi-colons to write-out
% these values):
sum(energy_sum)        % energy conservation total error
sum((mu_not*pi*Fs)+Fd) % total incident insolation (W m-2)
sum(sum(F_abs))        % total energy absorbed by all snow layers
sum(F_btm_net)         % total energy absorbed by underlying substrate


% Spectrally-integrated solar, visible, and NIR albedos:
alb     = sum(flx_slr.*albedo)./sum(flx_slr);

alb_vis = sum(flx_slr(1:vis_max_idx).*albedo(1:vis_max_idx))/...
          sum(flx_slr(1:vis_max_idx));

alb_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*albedo(vis_max_idx+1:nir_max_idx))/...
          sum(flx_slr(vis_max_idx+1:nir_max_idx));


% Diagnostic for comparing 470-band solutions with 5-band solutions:
%  Spectrally-integrated 5-band albedo:
if (BND_TYP == 1)
    bnd1a=1;
    bnd1b=40;
    bnd2a=41;
    bnd2b=70;
    bnd3a=71;
    bnd3b=90;
    bnd4a=91;
    bnd4b=120;
    bnd5a=121;
    bnd5b=470;
  
    bnd6a=1;
    bnd6b=40;
    bnd7a=41;
    bnd7b=89;
    bnd8a=90;
    bnd8b=470;
  
    alb1 = sum(flx_slr(bnd1a:bnd1b).*albedo(bnd1a:bnd1b))/sum(flx_slr(bnd1a:bnd1b));
    alb2 = sum(flx_slr(bnd2a:bnd2b).*albedo(bnd2a:bnd2b))/sum(flx_slr(bnd2a:bnd2b));
    alb3 = sum(flx_slr(bnd3a:bnd3b).*albedo(bnd3a:bnd3b))/sum(flx_slr(bnd3a:bnd3b));
    alb4 = sum(flx_slr(bnd4a:bnd4b).*albedo(bnd4a:bnd4b))/sum(flx_slr(bnd4a:bnd4b));
    alb5 = sum(flx_slr(bnd5a:bnd5b).*albedo(bnd5a:bnd5b))/sum(flx_slr(bnd5a:bnd5b));
    alb6 = sum(flx_slr(bnd6a:bnd6b).*albedo(bnd6a:bnd6b))/sum(flx_slr(bnd6a:bnd6b));
    alb7 = sum(flx_slr(bnd7a:bnd7b).*albedo(bnd7a:bnd7b))/sum(flx_slr(bnd7a:bnd7b));
    alb8 = sum(flx_slr(bnd8a:bnd8b).*albedo(bnd8a:bnd8b))/sum(flx_slr(bnd8a:bnd8b));
end;

% Spectrally-integrated VIS and NIR total snowpack absorption:
abs_vis = sum(flx_slr(1:vis_max_idx).*(1-albedo(1:vis_max_idx)));
abs_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*(1-albedo(vis_max_idx+1:nir_max_idx)));



% Output diagnostics:
%  1. The basics:
data_out(:,1) = wvl;            % spectral wavelength bands (um)
data_out(:,2) = albedo;         % spectral hemispheric albedo
data_out(1,3) = alb;            % solar broadband albedo
data_out(2,3) = alb_vis;        % visible (0.3-0.7um) albedo
data_out(3,3) = alb_nir;        % near-IR (0.7-5.0um) albedo
data_out(4,3) = sum(F_abs_slr); % total radiative absorption by all snow layers (not including underlying substrate)
data_out(5,3) = F_abs_slr(1);   % top layer solar absorption
if (nbr_lyr > 1)
    data_out(6,3) = F_abs_slr(2); % 2nd layer solar absorption
else
    data_out(6,3) = NaN;          % 2nd layer solar absorption
end;
data_out(:,5)  = sum(F_abs,2);  % spectral absorption

% more detail:
if (BND_TYP == 1)
    % different band-weighted albedos:
    data_out(1,4) = alb1;
    data_out(2,4) = alb2;
    data_out(3,4) = alb3;
    data_out(4,4) = alb4;
    data_out(5,4) = alb5;
    data_out(7,3) = alb6;
    data_out(8,3) = alb7;
    data_out(9,3) = alb8;
end;

data_out(6,4) = F_abs_slr(1); % top layer solar absorption
data_out(7,4) = F_abs_vis(1); % top layer VIS absorption
data_out(8,4) = F_abs_nir(1); % top layer NIR absorption

if (nbr_lyr > 1)
    data_out(9,4) = F_abs_slr(2);
    data_out(10,4) = F_abs_vis(2);
    data_out(11,4) = F_abs_nir(2);
    if (nbr_lyr > 2)
        data_out(12,4) = F_abs_slr(3);
        data_out(13,4) = F_abs_vis(3);
        data_out(14,4) = F_abs_nir(3);
        if (nbr_lyr > 3)
            data_out(15,4) = F_abs_slr(4);
            data_out(16,4) = F_abs_vis(4);
            data_out(17,4) = F_abs_nir(4);
            if (nbr_lyr > 4)
                data_out(18,4) = F_abs_slr(5);
                data_out(19,4) = F_abs_vis(5);
                data_out(20,4) = F_abs_nir(5);
            end;
        end;
    end;
end;

data_out(18,4) = F_abs_btm;     % solar absorption by underlying surface
data_out(19,4) = F_abs_vis_btm; % VIS absorption by underlying surface
data_out(20,4) = F_abs_nir_btm; % NIR absorption by underlying surface
data_out(21,4) = sum((mu_not*pi*Fs))+sum(Fd);  % total downwelling energy on upper boundary

% add by che for diagnose snow optical properties by internally mixed BC
data_out(:,6) = omega_snw(:,1); % spectral snow single-scattering albedo
data_out(:,7) = 1.0 - omega_snw(:,1); % snow single-scattering coalbedo
data_out(:,8) = g_snw(:,1); % spectral snow asymmetry factor

% plot modeled albedo:
if (1==0)
    plot(wvl,albedo,'k-');
    axis([0.3 2.5 0 1]);
    grid on;
end;

