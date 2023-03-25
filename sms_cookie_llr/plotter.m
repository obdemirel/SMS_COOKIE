clear all
rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);
ims = 2;


load('cine_r1_reg_paper_1.mat')
a1 = squeeze(mod_kspace(:,:,:,3,ims));
a1 = fliplr(rssq(ifft2(a1)));
a1 = a1./max(a1(:));

load('C:\Users\omer\Desktop\s10\prop_without\cine_r2_wo_reg_paper_1.mat')
a2 = squeeze(mod_kspace(:,:,:,3,ims));
a2 = fliplr(rssq(ifft2(a2)));
a2 = a2./max(a2(:));


figure, imshow(a1(46:112,42:98),[0 0.8])
figure, imshow(a2(46:112,42:98),[0 0.8])
figure, imshow(abs(a1(46:112,42:98)-a2(46:112,42:98)),[]), colorbar