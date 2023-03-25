function [ro_images] = sense_maker(recon,sense_maps,ims,slice_R)

ksb = size(recon,1)/slice_R;

[m,n,no_c,ims] = size(recon);

for ss=1:size(recon,4)
    concatenated_imags = encoder(recon(:,:,:,ss),m,n,no_c,slice_R,sense_maps,ksb);
    for slis = 1:slice_R
        ro_images(:,:,slis,ss) = concatenated_imags(ksb*(slis-1) + 1:slis*ksb,:);
    end
end

end