clear all
clc

%parpool(maxNumCompThreads)
%parpool(8)

for scans = 1:8;%1;%1:9
  

dirac = strcat('/home/range3-raid1/moeller-data/hcp_ret_7T_restore/TMP_DATA/20150822/SCAN',num2str(scans));

allFiles = dir([dirac '.mat']);


%files = dir(['/home/range3-raid1/moeller-data/hcp_ret_7T_restore/20150613/SCAN1'])
files = dir(dirac);
inner_counter = 1;
for slice_set = 5:21
    tic
load_name = files(slice_set).name;
load(strcat(dirac,'/',load_name))

load_dir = '/home/daedalus2-data2/icarus/Burak_Files/NORDIC_maps/test_20150822/';
load_name = strcat('sense_maps_',load_name(1:end-4),'_SCAN',num2str(scans));


load(strcat(load_dir,load_name))

%kspace_small = squeeze(DATA_mb_NORDIC(:,:,1,1,1,:,100));
%kspace = zeros(130,130,32,'single');
%kspace(:,18:4:end,:) = kspace_small;

kspace_small = squeeze(DATA_mb(:,:,1,1,1,:,98:101));
[mm,nn,no_c,ims] = size(kspace_small);
kspace_R2 = zeros(mm,130,no_c,ims,'single');
kspace_R2(:,18:2:end,:,:) = kspace_small;

[m,n,no_c] = size(kspace_R2);

sense_maps = single(sense_maps);

acs = OUT3.ACS;


%%%% SPSG MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parpool(4) %parpool(maxNumCompThreads)

[m,n,no_c,ims]= size(kspace_R2);
slice_R = 5;
PE_R = 2;
plotter = 0;

sg_shifter_tukey(m,n,acs)
disp('ACSs are Ready!')

tic
sg_kernel_main_sp(slice_R,PE_R,kspace_R2,'inital_ind')
toc
disp('SpSg Kernels are Ready!')

tic
sgrecon = sg_rec(slice_R,PE_R,ims,'sg_kernel_55',kspace_R2);
toc
disp('SpSg 1st Part is Ready!')

tic
finalrecon = sg_fin(slice_R,PE_R,ims,sgrecon,'inital_ind',kspace_R2);
toc
disp('SpSg 2nd Part is Ready!')


tic
[sense1_images] = sense_maker(finalrecon,sense_maps);
toc
disp('Sense1 Images are Ready!')

%%%%%%%%% end of SPSG MAIn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sense1_images_all(:,:,:,:,inner_counter) = sense1_images;

inner_counter = inner_counter+1
toc
end

save_dir = '/home/naxos2-raid1/omer/spsg_to_server/results/';
save_name = strcat('subject_',num2str(scans));

%save(fullfile(save_dir,save_name),'sense1_images_all','-v7.3')

save(save_name,'sense1_images_all','-v7.3')
end

delete(gcp)