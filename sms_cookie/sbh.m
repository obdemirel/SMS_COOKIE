close all
clear all
clc

recon_align

rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);


figure, subplot(1,5,1), imshow(flipud(rssq(ifftshift(ifft2(squeeze(ACS_sl(:,:,:,1))),1))),[0 1e-5])
subplot(1,5,2), imshow(flipud(rssq(ifftshift(ifft2(squeeze(ACS_sl(:,:,:,2))),1))),[0 1e-5])
subplot(1,5,3), imshow(flipud(rssq(ifftshift(ifft2(squeeze(ACS_sl(:,:,:,3))),1))),[0 1e-5])
subplot(1,5,4), imshow(flipud(rssq(ifftshift(ifft2(squeeze(ACS_sl(:,:,:,4))),1))),[0 1e-5])
subplot(1,5,5), imshow(flipud(rssq(ifftshift(ifft2(squeeze(ACS_sl(:,:,:,5))),1))),[0 1e-5])


[m,n,no_c,ims] = size(kspace);
[mm,nn,no_c,no_slices] = size(ACS_sl);

%% kspace mid is at 118-89 shift accordingly



%% 1+4+7
ACS1_kspace = zeros(m,n,no_c,'single');
ACS1_kspace(:,76-34+1:76+30,:) = ACS_sl(:,:,:,1);
ACS1 = circshift(ACS1_kspace,[0 0 0]);

phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';
phas2 = repmat(phases2,[1 n no_c]);
ACS11 = ACS1 .* phas2;
ACS1_final = ACS11(:,76-34+1:76+30,:);
ACS11_final = circshift(ACS1_final,[0 0 0]);
c = sum(abs((ACS11_final(:,:,:))),3);
mc = max(c(:));
[c1,c2] = find(c==mc);

% %%% for 2
ACS1_kspace = zeros(m,n,no_c,'single');
ACS1_kspace(:,76-34+1:76+30,:) = ACS_sl(:,:,:,2);
ACS1 = circshift(ACS1_kspace,[0 0 0]);

phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';
phas2 = repmat(phases2,[1 n no_c]);
ACS11 = ACS1 .* phas2;
ACS1_final = ACS11(:,76-34+1:76+30,:);
ACS2_final = circshift(ACS1_final,[0 0 0]);
c = sum(abs((ACS2_final(:,:,:))),3);
mc = max(c(:));
[c1,c2] = find(c==mc);

% %%% for 3
ACS1_kspace = zeros(m,n,no_c,'single');
ACS1_kspace(:,76-34+1:76+30,:) = ACS_sl(:,:,:,3);
ACS1 = circshift(ACS1_kspace,[0 0 0]);

phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';
phas2 = repmat(phases2,[1 n no_c]);
ACS11 = ACS1 .* phas2;
ACS1_final = ACS11(:,76-34+1:76+30,:);
ACS3_final = circshift(ACS1_final,[0 0 0]);
c = sum(abs((ACS3_final(:,:,:))),3);
mc = max(c(:));
[c1,c2] = find(c==mc);

% %%% for 4
ACS1_kspace = zeros(m,n,no_c,'single');
ACS1_kspace(:,76-34+1:76+30,:) = ACS_sl(:,:,:,4);
ACS1 = circshift(ACS1_kspace,[0 0 0]);

phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';
phas2 = repmat(phases2,[1 n no_c]);
ACS11 = ACS1 .* phas2;
ACS1_final = ACS11(:,76-34+1:76+30,:);
ACS4_final = circshift(ACS1_final,[0 0 0]);
c = sum(abs((ACS4_final(:,:,:))),3);
mc = max(c(:));
[c1,c2] = find(c==mc);

% %%% for 5
ACS1_kspace = zeros(m,n,no_c,'single');
ACS1_kspace(:,76-34+1:76+30,:) = ACS_sl(:,:,:,5);
ACS1 = circshift(ACS1_kspace,[0 0 0]);

phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';
phas2 = repmat(phases2,[1 n no_c]);
ACS11 = ACS1 .* phas2;
ACS1_final = ACS11(:,76-34+1:76+30,:);
ACS5_final = circshift(ACS1_final,[0 0 0]);
c = sum(abs((ACS5_final(:,:,:))),3);
mc = max(c(:));
[c1,c2] = find(c==mc);


save_name = 'initial_cine_acs';
ACS1_final = ACS11_final;
save(save_name,'ACS1_final','ACS2_final','ACS3_final','ACS4_final','ACS5_final')

figure, subplot(1,5,1), imshow(flipud(rssq((ifft2(squeeze(ACS1_final(:,:,:)))))),[0 1e-5])
subplot(1,5,2), imshow(flipud(rssq((ifft2(squeeze(ACS2_final(:,:,:)))))),[0 1e-5])
subplot(1,5,3), imshow(flipud(rssq((ifft2(squeeze(ACS3_final(:,:,:)))))),[0 1e-5])
subplot(1,5,4), imshow(flipud(rssq((ifft2(squeeze(ACS4_final(:,:,:)))))),[0 1e-5])
subplot(1,5,5), imshow(flipud(rssq((ifft2(squeeze(ACS5_final(:,:,:)))))),[0 1e-5])
% 
%% DO the same for the MB Data
nkspace = zeros(m,n,no_c,ims,'single');
% mid1 = circshift(squeeze(kspace(:,:,:,selected_slab,:)),[0 6 0 0]);
nkspace = kspace;% .* phas;
nkspace = circshift(nkspace, [0 0 0]);
phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).'; 
phas2 = repmat(phases2,[1 n no_c ims]);
mnkspace = nkspace .* phas2;
new_kspace = circshift(mnkspace,[0 0 0]);

c = sum(sum(abs(new_kspace),4),3);
mc = max(c(:));
[c1,c2] = find(c==mc);

% % 
save_name = 'initial_cine_kspace';
save(save_name,'new_kspace')
% 
