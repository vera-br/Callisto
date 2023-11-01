% SPICEcookbook.m% Aug 17, 2010% stefano.campagnola@missionanalysis.org
% Read SPICEcookbook.pdf'
% -----------------------------------------------------------

% This is the folder with the kernels
% () Replace this string with your [datafolder] path!
datafolder = 'JUICE/JUICE/kernels/';

addpath('mice/mice/lib');
addpath('mice/mice/src/mice')

% Load leapsecond kernel
KERNELS_TO_LOAD = ["ck/juice_sc_crema_5_1_150lb_23_1_default_v01.bc"
"ck/juice_sc_crema_5_1_150lb_23_1_comms_v01.bc"
"ck/juice_sc_crema_5_1_150lb_23_1_conjctn_v01.bc"
"ck/juice_sc_crema_5_1_150lb_23_1_flybys_v01.bc"
"ck/juice_sc_crema_5_1_150lb_23_1_baseline_v02.bc"
"ck/juice_lpbooms_f160326_v01.bc"
"ck/juice_magboom_f160326_v04.bc"
"ck/juice_majis_scan_zero_v02.bc"
"ck/juice_swi_scan_zero_v02.bc"
"ck/juice_sa_crema_5_1_150lb_23_1_default_v01.bc"
"ck/juice_sa_crema_5_1_150lb_23_1_baseline_v03.bc"
"ck/juice_mga_crema_5_1_150lb_23_1_default_v01.bc"
"ck/juice_mga_crema_5_1_150lb_23_1_baseline_v03.bc"
"fk/juice_v36.tf"
"fk/juice_sci_v17.tf"
"fk/juice_ops_v08.tf"
"fk/juice_dsk_surfaces_v09.tf"
"fk/juice_roi_v02.tf"
"fk/juice_events_crema_5_1_150lb_23_1_v02.tf"
"fk/rssd0002.tf"
"fk/earth_topo_050714.tf"
"fk/earthfixediau.tf"
"fk/estrack_v04.tf"
"dsk/juice_europa_plasma_torus_v03.bds"
"dsk/juice_io_plasma_torus_v05.bds"
"dsk/juice_jup_ama_gos_ring_v02.bds"
"dsk/juice_jup_halo_ring_v04.bds"
"dsk/juice_jup_main_ring_v04.bds"
"dsk/juice_jup_the_ring_ext_v01.bds"
"dsk/juice_jup_the_gos_ring_v02.bds"
"dsk/juice_sc_bus_v06.bds"
"dsk/juice_sc_gala_v02.bds"
"dsk/juice_sc_janus_v02.bds"
"dsk/juice_sc_jmc1_v02.bds"
"dsk/juice_sc_jmc2_v02.bds"
"dsk/juice_sc_lpb1_v02.bds"
"dsk/juice_sc_lpb2_v02.bds"
"dsk/juice_sc_lpb3_v02.bds"
"dsk/juice_sc_lpb4_v02.bds"
"dsk/juice_sc_rwi_v02.bds"
"dsk/juice_sc_scm_v02.bds"
"dsk/juice_sc_mag_v06.bds"
"dsk/juice_sc_majis_v02.bds"
"dsk/juice_sc_mga_apm_v04.bds"
"dsk/juice_sc_mga_dish_v04.bds"
"dsk/juice_sc_pep_jdc_v02.bds"
"dsk/juice_sc_pep_jei_v02.bds"
"dsk/juice_sc_pep_jeni_v02.bds"
"dsk/juice_sc_pep_jna_v02.bds"
"dsk/juice_sc_pep_nim_v02.bds"
"dsk/juice_sc_rimemx_v02.bds"
"dsk/juice_sc_rimepx_v02.bds"
"dsk/juice_sc_sapy_v01.bds"
"dsk/juice_sc_samy_v01.bds"
"dsk/juice_sc_str1_v02.bds"
"dsk/juice_sc_str2_v02.bds"
"dsk/juice_sc_str3_v02.bds"
"dsk/juice_sc_swi_v02.bds"
"ik/juice_gala_v05.ti"
"ik/juice_janus_v08.ti"
"ik/juice_jmc_v02.ti"
"ik/juice_jmag_v02.ti"
"ik/juice_majis_v08.ti"
"ik/juice_navcam_v01.ti"
"ik/juice_pep_v13.ti"
"ik/juice_radem_v02.ti"
"ik/juice_rime_v04.ti"
"ik/juice_rpwi_v03.ti"
"ik/juice_str_v01.ti"
"ik/juice_swi_v07.ti"
"ik/juice_uvs_v06.ti"
"ik/juice_aux_v02.ti"
"lsk/naif0012.tls"
"pck/pck00011.tpc"
"pck/de-403-masses.tpc"
"pck/gm_de431.tpc"
"pck/inpop19a_moon_pa_v01.bpc"
"pck/earth_070425_370426_predict.bpc"
"pck/juice_jup010.tpc"
"pck/juice_roi_v01.tpc"
"sclk/juice_fict_160326_v02.tsc"
"spk/juice_sci_v04.bsp"
"spk/juice_struct_v19.bsp"
"spk/juice_cog_v00.bsp"
"spk/juice_roi_v02.bsp"
"spk/mar085_20200101_20400101.bsp"
"spk/earthstns_fx_050714.bsp"
"spk/estrack_v04.bsp"
"spk/jup365_19900101_20500101.bsp"
"spk/jup343_19900101_20500101.bsp"
"spk/jup344-s2003_j24_19900101_20500101.bsp"
"spk/jup346_19900101_20500101.bsp"
"spk/de432s.bsp"
"spk/noe-5-2017-gal-a-reduced_20200101_20380902.bsp"
"spk/juice_crema_5_1_150lb_23_1_v01.bsp"];

for i = 1:length(KERNELS_TO_LOAD)

    chars_KERNELS = convertStringsToChars(KERNELS_TO_LOAD(i));
    cspice_furnsh([datafolder, chars_KERNELS]);

end

% Retrieve some constants:
% () For more info, read the documentation on cspice_bodvcd
kJ = cspice_bodvcd(5, 'GM', 10); % GM of Jupiter with 10 significant
% digits
rJ = cspice_bodvcd(503, 'RADII', 10);% Radius of Ganymede
% Choose a date and convert it into ephemeris time
% () For more info, read the documentation on cspice_str3et
hour = 3600;
date0 = '24 Jun 2034 00:45:00 UTC';
et0 = cspice_str2et(date0);
et_R = et0:3600/60:(et0 + 9 * hour);

% Ephemerides of Callisto w.r.t. Jupiter barycenter
% on the ecliptic J2000 at the times et_R
% () For more info, read the documentation on cspice_spkezr
callisto_data = cspice_spkezr('504', et_R, 'JUICE_CALLISTO_PHI_ORB', 'NONE', '599');
JUICE_pos = cspice_spkezr('-28', et_R, 'JUICE_CALLISTO_PHI_ORB', 'NONE', '504');
juice_jupiter = cspice_spkezr('-28', et_R, 'JUICE_CALLISTO_PHI_ORB', 'NONE', '599');
juice_jupiter_mag = cspice_spkezr('-28', et_R, 'JUPITER_MAG_VIP4', 'NONE', '599');
callisto_jupiter_mag = cspice_spkezr('504', et_R, 'JUPITER_MAG_VIP4', 'NONE', '599');


% Convert an ephemeris time back into calendar format% () For more info, read the documentation on cspice_et2utc
et1 = et0+7*86400;
date1 = cspice_et2utc(et1,'C',6);

writematrix([callisto_data;et_R], 'spice_data/callisto_wrt_jupiter_C21.csv');
writematrix([JUICE_pos;et_R], 'spice_data/juice_wrt_callisto_C21.csv');
writematrix([juice_jupiter;et_R], 'spice_data/juice_wrt_jupiter_C21.csv');
writematrix([juice_jupiter_mag;et_R], 'spice_data/juice_wrt_jupiter_mag_C21.csv');
writematrix([callisto_jupiter_mag;et_R], 'spice_data/callisto_wrt_jupiter_mag_C21.csv');

% \C’ = ’calendar’, with 6 digits
% Clear the memory
cspice_kclear