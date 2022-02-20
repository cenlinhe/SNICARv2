% Driver routine for SNICAR.  Also see commenting in snicar8d_v2.m
%
% Original version: written by Mark Flanner (flanner@umich.edu) Flanner et al. 2007 JGR
% Current version: Updated by Cenlin He (cenlinhe@ucar.edu) He et al. 2018 ACP
% Original version assumes aerosol externally mixed with spherical snow grains 
% Current version includes nonspherical snow grain & BC-snow internal mixing 
% based on parameterizations from [He et al. 2017,J. Climate, doi:10.1175/JCLI-D-17-0300.1]
% Dust-snow internal mixing parametriztaions from (He et al., 2019 JAMES)
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
% mss_cnc_dstint: mass mixing ratio of dust species (internally mixed with snow)
%                 (units of parts per billion, ng g-1) add by che
% mss_cnc_dstext: mass mixing ratio of dust species (externally mixed with snow)
%                 (units of parts per billion, ng g-1) add by che
% fl_sot1:      name of file containing optical properties for BC species 1
% fl_sot2:      name of file containing optical properties for BC species 2
% fl_dst1:      name of file containing optical properties for dust species 1
% fl_dst2:      name of file containing optical properties for dust species 2
% fl_dst3:      name of file containing optical properties for dust species 3
% fl_dst4:      name of file containing optical properties for dust species 4
% fl_ash1:      name of file containing optical properties for ash species 1
%
%
%%%%% che %%%%%% Options added by Cenlin He, Apr. 2018
% SNO_SHP:      1=sphere (SNICAR default shape)
%               2=spheroid; 3=hexagonal plate; 4=koch snowflake;
% SNO_fs:       Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
%               0=use recommended default value (He et al. 2017 JC);
%               others(0<fs<1)= use user-specified value 
%               only activated when SNO_SHP > 1 (i.e. nonspherical)
% SNO_AR:       Aspect ratio: ratio of grain width to length
%               0=use recommended default value (He et al. 2017 JC);
%               others(0.1<fs<20)= use user-specified value
%               only activated when SNO_SHP > 1 (i.e. nonspherical)
% BC_SNO_MIX:   1= only BC-snow external mixing (SNICAR default)
%               2= BC-snow internal mixing (He et al. 2017 JC)
% DST_SNO_MIX:  1= only dust-snow external mixing (SNICAR default)
%               2= dust-snow internal mixing (He et al. 2019 JAMES)
%%%%% che %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;



% RADIATIVE TRANSFER CONFIGURATION:
BND_TYP  = 1;        % 1= 470 spectral bands
DIRECT   = 0;        % 1= Direct-beam incident flux, 0= Diffuse incident flux
APRX_TYP = 3;        % 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
DELTA    = 1;        % 1= Apply Delta approximation, 0= No delta
                     
% COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM
coszen   = 1;

% REFLECTANCE OF SURFACE UNDERLYING SNOW:
%   Value is applied to all wavelengths.
%   User can also specify spectrally-dependent ground albedo
%   internally in snicar8d.m
R_sfc    = 0;

% SNOW LAYER THICKNESSES (array) (units: meters):
dz       = [1000]; 
 
nbr_lyr  = length(dz);  % number of snow layers

% SNOW DENSITY OF EACH LAYER (units: kg/m3)
rho_snw(1:nbr_lyr) = 150;  

% SNOW EFFECTIVE GRAIN SIZE FOR EACH LAYER (units: microns):
rds_snw(1:nbr_lyr) = 100;


%%%% che %%%%%
% Options added by Cenlin He for nonspherical snow & BC-snow internal
% mixing based on He et al. (2017) parameterization 
% He et al. 2017: J. Climate, 30, doi:10.1175/JCLI-D-17-0300.1

SNO_SHP(1:nbr_lyr)  = 1;    % Snow grain shape option
                            % 1=sphere (SNICAR default shape)
                            % 2=spheroid; 3=hexagonal plate; 4=koch snowflake;
                            
SNO_fs(1:nbr_lyr)   = 0;    % Shape factor: ratio of nonspherical grain effective radii to that of sphere with the same volume.
                            % 0=use recommended default value (He et al. 2017 JC);
                            % others(0<fs<1)= use user-specified value
                            % only activated when SNO_SHP > 1 (i.e. nonspherical)
                            
SNO_AR(1:nbr_lyr)   = 0;    % Aspect ratio: ratio of grain width to length
                            % 0=use recommended default value (He et al. 2017 JC);
                            % others(0.1<fs<20)= use user-specified value
                            % only activated when SNO_SHP > 1 (i.e. nonspherical)
                            
BC_SNO_MIX(1:nbr_lyr) = 1;  % BC-snow mixing option (He et al., 2017 JC)
                            % 1= only BC-snow external mixing (SNICAR default)
                            % 2= BC-snow internal mixing

DST_SNO_MIX(1:nbr_lyr) = 1; % dust-snow mixing option (He et al., 2019 JAMES)
                            % 1= only dust-snow external mixing (SNICAR default)
                            % 2= dust-snow internal mixing
%%%% che %%%%%

  
% NUMBER OF AEROSOL SPECIES IN SNOW (ICE EXCLUDED)
%  Species numbers (used in snicar8d.m) are:
%    1: uncoated black carbon
%    2: coated black carbon
%    3: dust size 1
%    4: dust size 2
%    5: dust size 3
%    6: dust size 4
%    7: volcanic ash
nbr_aer = 7;

% PARTICLE MASS MIXING RATIOS (units: ng(species)/g(ice), or ppb)
mss_cnc_dstint(1:nbr_lyr)= 0.0;    % che: dust internally mixed with snow
mss_cnc_dstext(1:nbr_lyr)= 0.0;    % che: dust externally mixed with snow
mss_cnc_sot3(1:nbr_lyr)  = 0.0;    % che: black carbon internally mixed with snow
mss_cnc_sot1(1:nbr_lyr)  = 0.0;    % uncoated black carbon
mss_cnc_sot2(1:nbr_lyr)  = 0.0;    % coated black carbon
mss_cnc_dst1(1:nbr_lyr)  = 0.0;    % dust species 1
mss_cnc_dst2(1:nbr_lyr)  = 0.0;    % dust species 2
mss_cnc_dst3(1:nbr_lyr)  = 0.0;    % dust species 3
mss_cnc_dst4(1:nbr_lyr)  = 0.0;    % dust species 4
mss_cnc_ash1(1:nbr_lyr)  = 0.0;    % volcanic ash species 1
 

 
% FILE NAMES CONTAINING MIE PARAMETERS FOR ALL AEROSOL SPECIES:
fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';
fl_dst1  = 'aer_dst_bln_20060904_01.nc';
fl_dst2  = 'aer_dst_bln_20060904_02.nc';
fl_dst3  = 'aer_dst_bln_20060904_03.nc';
fl_dst4  = 'aer_dst_bln_20060904_04.nc';
fl_ash1  = 'volc_ash_mtsthelens_20081011.nc';



% call SNICAR with these inputs:
% Original SNICAR by calling "snicar8d"
% che: Updated version by calling "snicar8d_v2"
% che: Updated version by calling "snicar8d_v2_dustint"
data_in = snicar8d_v2_dustint(BND_TYP, DIRECT, APRX_TYP, DELTA, coszen, R_sfc, ...
                   dz, rho_snw, rds_snw, nbr_aer, mss_cnc_sot1, ...
                   mss_cnc_sot2, mss_cnc_dst1, mss_cnc_dst2, ...
                   mss_cnc_dst3, mss_cnc_dst4, mss_cnc_ash1, fl_sot1, ...
                   fl_sot2, fl_dst1, fl_dst2, fl_dst3, fl_dst4, fl_ash1, ...
                   SNO_SHP, BC_SNO_MIX, SNO_fs, SNO_AR, mss_cnc_sot3, ...
                   DST_SNO_MIX, mss_cnc_dstint, mss_cnc_dstext); 


% process input data:
% (see description of data_out at the end of snicar8d_v2.m for more info)
wvl         = data_in(:,1);   % wavelength grid
albedo      = data_in(:,2);   % spectral albedo
alb_slr     = data_in(1,3);   % broadband albedo (0.3-5.0um)
alb_vis     = data_in(2,3);   % visible albedo (0.3-0.7um)
alb_nir     = data_in(3,3);   % near-IR albedo (0.7-5.0um)
flx_abs_snw = data_in(4,3);   % total radiative absorption by all snow layers (not including underlying substrate)

flx_abs(1)     = data_in(6,4); % top layer solar absorption
flx_vis_abs(1) = data_in(7,4); % top layer VIS absorption
flx_nir_abs(1) = data_in(8,4); % top layer NIR absorption

% make a plot of spectrally-resolved albedo:
plot(wvl,albedo,'linewidth',3);
xlabel('Wavelength (\mum)','fontsize',20);
ylabel('Albedo','fontsize',20);
xlim([0.3,5]);
ylim([0,1]);
set(gca,'xscale','log');
set(gca,'xtick',0.3:0.5:5.0,'fontsize',16);
set(gca,'ytick',0:0.1:1.0,'fontsize',16);
grid on;


