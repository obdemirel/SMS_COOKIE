clear all

load data
dummy = kspace_R2;
kspace_R2(:,:,:,2) = dummy;
kspace_R2(:,:,:,3) = dummy;
kspace_R2(:,:,:,4) = dummy;

%parpool(4) %parpool(maxNumCompThreads)

[m,n,no_c,ims]= size(kspace_R2);
slice_R = 5;
PE_R = 2;
plotter = 1;

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

load sense_maps

tic
[sense1_images] = sense_maker(finalrecon,sense_maps);
toc
disp('Sense1 Images are Ready!')


if(plotter==1)
im1 = abs(sense1_images(:,:,1,1));
im2 = abs(sense1_images(:,:,2,1));
im3 = abs(sense1_images(:,:,3,1));
im4 = abs(sense1_images(:,:,4,1));
im5 = abs(sense1_images(:,:,5,1));
im1 = im1./max(im1(:));
im2 = im2./max(im2(:));
im3 = im3./max(im3(:));
im4 = im4./max(im4(:));
im5 = im5./max(im5(:));

figure, imshow(cat(2,im1,im2,im3,im4,im5),[]), ylabel('SMS5 x R2 SpSg Recon')
end
