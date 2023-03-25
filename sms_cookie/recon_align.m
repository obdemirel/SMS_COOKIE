clear all
%path(path, '/home/daedalus1-raid1/akcakaya-group-data/ScannerData/Data/Volunteer/2017-08-24 Ac_CMR_Vol_022 - Contrast/MBRS_cine');

% load PatRef_Cine_630
load PatRef_Cine_632
[m,n,no_c,no_MB] = size(Coils_kspace);
Coils_kspace = Coils_kspace(:,:,:,[1 4 2 5 3]);
Coils_kspace = cat(2, Coils_kspace, zeros(m,1,no_c,no_MB,'single'));
[m,n,no_c,no_MB] = size(Coils_kspace);


for ind = 1:no_MB
    phases = exp(sqrt(-1)*2*pi*ind/no_MB*(0:n-1));
    phas = repmat(phases,[m 1 no_c]);
    ACS_sl(:,:,:,ind) = Coils_kspace(:,:,:,no_MB - ind+1).* phas;
end


m_acs = m;
n_acs = n;

load meas_MID00633_FID20455_FLASH_CINE_MBR5_8int_8seg_Retro
% load meas_MID00636_FID20458_FLASH_CINE_MBR5_noRS_6mm_Retro

% shift the center correctly
[m,n,no_c,no_ims] = size(kSpaceSorted);
kspace = permute(kSpaceSorted, [1 2 4 3]);
kspace = kspace(:,:,:,2:end-1);

