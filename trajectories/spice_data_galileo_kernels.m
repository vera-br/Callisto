% SPICEcookbook.m% Aug 17, 2010% stefano.campagnola@missionanalysis.org
% Read SPICEcookbook.pdf'
% -----------------------------------------------------------

% This is the folder with the kernels
% () Replace this string with your [datafolder] path!
datafolder = './GLL/kernels/';

addpath('mice/lib');
addpath('mice/src/mice')

% Load leapsecond kernel
% ck - fk - dsk -ik - lsk - pck - sclk - spk
KERNELS_TO_LOAD = ["ck/gll_plt_pre_1990_v00.bc"
"ck/gll_plt_pre_1991_v00.bc"
"ck/gll_plt_pre_1992_v00.bc"
"ck/gll_plt_pre_1993_v00.bc"
"ck/gll_plt_pre_1994_v00.bc" 
"ck/gll_plt_pre_1995_v00.bc" 
"ck/gll_plt_pre_1996_v00.bc"
"ck/gll_plt_pre_1997_v00.bc" 
"ck/gll_plt_pre_1998_v00.bc"
"ck/gll_plt_pre_1999_v00.bc" 
"ck/gll_plt_pre_2000_v00.bc" 
"ck/gll_plt_pre_2001_v00.bc"
"ck/gll_plt_rec_1990_tav_v00.bc" 
"ck/gll_plt_rec_1991_tav_v00.bc" 
"ck/gll_plt_rec_1993_tav_v00.bc"
"ck/gll_plt_rec_1994_tav_v00.bc" 
"ck/gll_plt_rec_1995_mav_v00.bc" 
"ck/gll_plt_rec_1995_tav_v00.bc" 
"ck/gll_plt_rec_1996_mav_v00.bc"
"ck/gll_plt_rec_1997_mav_v00.bc" 
"ck/gll_plt_rec_1998_mav_v00.bc" 
"ck/gll_plt_rec_1999_mav_v00.bc"
"ck/gll_plt_rec_2000_mav_v00.bc"
"ck/gll_plt_rec_2001_mav_v00.bc"
"ck/gll_plt_rec_2002_mav_v00.bc"
"ik/gll24001.ti" 
"ik/gll25001.ti" 
"ik/gll28001.ti" 
"ik/gll29001.ti" 
"ik/gll34001.ti" 
"ik/gll36001.ti"
"ik/gll37001.ti"
"lsk/mk98264a.tls"
"lsk/naif0008.tls"
"pck/mips_010314.tpc" 
"pck/pck00007.tpc"
"pck/pk96030a.tpc"
"sclk/mk00062a.tsc"
"spk/gll_951120_021126_raj2007.bsp"
"spk/gll_951120_021126_raj2007.bsp.lbl" 
"spk/gll_951120_021126_raj2021.bsp"
"spk/gll_951120_021126_raj2021.bsp.lbl" 
"spk/s000131a.bsp"
"spk/s000131a.bsp.lbl"
"spk/s030916a.bsp"
"spk/s030916a.bsp.lbl"
"spk/s960730a.bsp"
"spk/s960730a.bsp.lbl"
"spk/s970311a.bsp"
"spk/s970311a.bsp.lbl" 
"spk/s980326a.bsp"
"spk/s980326a.bsp.lbl"];

KERNELS_TO_LOAD3 = ["ck/juice_sc_crema_5_1_150lb_23_1_default_v01.bc"
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

for i = 1:length(KERNELS_TO_LOAD3)

    chars_KERNELS = convertStringsToChars(KERNELS_TO_LOAD3(i));
    cspice_furnsh(['./JUICE/kernels/', chars_KERNELS]);

end

closest_approach_times = ["1996-11-04 13:34:28 UTC"
"1997-06-25 13:47:50.279000 UTC"];

% Retrieve some constants:
% () For more info, read the documentation on cspice_bodvcd
kJ = cspice_bodvcd(5, 'GM', 10); % GM of Jupiter with 10 significant
% digits
rJ = cspice_bodvcd(503, 'RADII', 10);% Radius of Ganymede

hour = 3600;
day = 86400;

for i = 1:length(closest_approach_times)

    date_CA = closest_approach_times(i);
    date_CA_i = convertStringsToChars(date_CA);
    et_CA = cspice_str2et(date_CA_i);
    et_0 = et_CA - 0.5 * day;
    et_R = et_0:3600/60:(et_0 + 1 * day);
    

    callisto_jupiter_SIII = cspice_spkezr('GALILEO ORBITER', et_R, 'JUICE_JUPITER_MAG_S3RH2009', 'NONE', '599');
    writematrix([callisto_jupiter_SIII;et_R], append('spice_data/galileo_wrt_jupiter_SIII_GK',string(i),'.csv'));

    galileo_callisto_cphio = cspice_spkezr('GALILEO ORBITER', et_R, 'JUICE_JUPITER_MAG_S3RH2009', 'NONE', '504');
    writematrix([galileo_callisto_cphio;et_R], append('spice_data/galileo_wrt_callisto_cphio_GK',string(i),'.csv'));
    


end
    
cspice_kclear
    
    
