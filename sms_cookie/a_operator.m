function y = a_operator(x, m, n, no_c)

x = reshape(x,m,n,no_c);

k1 = x(1:164,:,:);
k2 = x(165:2*164,:,:);
k3 = x(329:3*164,:,:);

% for co =1:no_c
% im1(:,:,co) = ifft2(k1(:,:,co));
% im2(:,:,co) = ifft2(k2(:,:,co));
% im3(:,:,co) = ifft2(k3(:,:,co));
% end

im1 = ifft2(k1);
im2 = ifft2(k2);
im3 = ifft2(k3);

im_all = [im1;im2;im3];

% for co =1:no_c
% y(:,:,co) = fft2(im_all(:,:,co));
% end

y = fft2(im_all);

y = y(:);

end