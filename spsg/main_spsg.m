%%%% SPSG MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parpool(4) %parpool(maxNumCompThreads)

[m,n,no_c,ims]= size(kspace);


tic
sg_kernel_main_sp(slice_R,PE_R,kspace,acs)
[data_ak_ind,kernel_r,kernel_s] = sg_kernel_main_sp(slice_R,PE_R,kspace,acs,[7,7])
toc
disp('SpSg Kernels are Ready!')

tic
sgrecon = sg_rec(slice_R,PE_R,ims,'sg_kernel_55',kspace);
toc
disp('SpSg 1st Part is Ready!')

tic
finalrecon = sg_fin(slice_R,PE_R,ims,sgrecon,'inital_ind',kspace);
toc
disp('SpSg 2nd Part is Ready!')


tic
[sense1_images] = sense_maker(finalrecon,sense_maps);
toc
disp('Sense1 Images are Ready!')


