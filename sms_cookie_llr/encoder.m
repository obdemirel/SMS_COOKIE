function [out] = encoder(x,m,n,no_c,slice_n,Coil_sensitivites,ksb)

kspace_to_im = @(x) ifft2c(x);
z_im = reshape(x,m,n,no_c);

for slice_no =1:slice_n
    slice_im = kspace_to_im(z_im((slice_no-1)*ksb + 1:slice_no*ksb,:,:));
    sense_images(:,:,slice_no) = sum(conj(Coil_sensitivites(:,:,:,slice_no)).*slice_im,3);
end

out =[];
for slis = 1:slice_n
    out = [out;sense_images(:,:,slis)];%;sense_images(:,:,4);sense_images(:,:,5)];
end
end

