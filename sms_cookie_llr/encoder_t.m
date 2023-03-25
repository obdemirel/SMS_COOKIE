function [out] = encoder_t(x,m,n,no_c,slice_n,Coil_sensitivites,ksb)

im_to_kspace = @(x) fft2c(x);
z_im1 = reshape(x,m,n);

for slis = 1:slice_n
    z_im(:,:,slis) = z_im1(ksb*(slis-1) + 1:slis*ksb,:);
end



for slice_no =1:slice_n
    new_z((slice_no-1)*ksb + 1:slice_no*ksb,:,:) = repmat(z_im(:,:,slice_no),[1 1 no_c]).*Coil_sensitivites(:,:,:,slice_no);
    z_kspace((slice_no-1)*ksb + 1:slice_no*ksb,:,:) = im_to_kspace(new_z((slice_no-1)*ksb + 1:slice_no*ksb,:,:));
end

out = [];
for slis = 1:slice_n
    out = [out;z_kspace(ksb*(slis-1) + 1:slis*ksb,:,:)];%;z_kspace(ksb*3 +1 :4*ksb,:,:);z_kspace(ksb*4 +1 :5*ksb,:,:)];
end
end

