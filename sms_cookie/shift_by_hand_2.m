clear all
% % 
%  load 10_MB_166_168_24ACS
% load MB_030_031;
% load 13_MB_089_094
load 16_MB_141_151
%kspace(:,2:2:end,:,:) = 0;s
 %load MB_166_168
% load 14_MB_030_031_24ACS
%% First apply the CAIPIRINHA shifts 
pre_processing_for_MB_data 
%%if you don't want the CAIPIRINHA  please uncomment the following 3 lines
% ACS_1 = squeeze(Coils_kspace(:,:,:,1));
% ACS_2 = squeeze(Coils_kspace(:,:,:,2));
% ACS_3 = squeeze(Coils_kspace(:,:,:,3));

slice_R = 3;
PE_R = 2;

rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);

[m,n,no_c,ims] = size(kspace);
[mm,nn,no_c] = size(ACS_1);

ACS1_kspace = zeros(m,n,no_c,'single');
ACS1_kspace(34:34+mm-1,45:45+nn-1,:) = ACS_1;
%checking the center position is correct or not
c = sum(abs((ACS1_kspace(:,:,:))),3);
mc = max(c(:));
[c1,c2] = find(c==mc);
ACS2_kspace = zeros(m,n,no_c,'single');
ACS2_kspace(34:34+mm-1,45:45+nn-1,:) = ACS_2;
% c = sum(abs(squeeze(ACS2_kspace(:,:,:))),3);
% mc = max(c(:));
% [c1,c2] = find(c==mc);
ACS3_kspace = zeros(m,n,no_c,'single');
ACS3_kspace(34:34+mm-1,45:45+nn-1,:) = ACS_3;
% c = sum(abs(squeeze(ACS3_kspace(:,:,:))),3);
% mc = max(c(:));
% [c1,c2] = find(c==mc);

phases = exp(sqrt(-1)*pi*(0:n-1));
phas = repmat(phases,[m 1 no_c]);
ACS1 = ACS1_kspace;% .* phas;
ACS2 = ACS2_kspace;% .* phas;
ACS3 = ACS3_kspace;% .* phas;
phases2 = (exp(sqrt(-1)*(pi+(2*pi/m))*(0:m-1))).';
%phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';%% intentionally add 2*pi/m for 1 pixel shift
%% if you want you can uncomment the line below and see what happens
% phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).';
phas2 = repmat(phases2,[1 n no_c]);
ACS11 = ACS1 .* phas2;
ACS22 = ACS2 .* phas2;
ACS33 = ACS3 .* phas2;

%% restore the zero padding
ACS1_final = ACS11(34:34+mm-1,45:45+nn-1,:); 
ACS2_final = ACS22(34:34+mm-1,45:45+nn-1,:);
ACS3_final = ACS33(34:34+mm-1,45:45+nn-1,:);

acs_im1 = ifft2(ACS1_final);
acs_im2 = ifft2(ACS2_final);
acs_im3 = ifft2(ACS3_final);

figure, imshow(rssq(acs_im1),[]), title('slice 1 ACS Image')
figure, imshow(rssq(acs_im2),[]), title('slice 2 ACS Image')
figure, imshow(rssq(acs_im3),[]), title('slice 3 ACS Image')

acs_im = [acs_im1;acs_im2;acs_im3];
%acs_im = [acs_im1(1:63,:,:);acs_im2;acs_im3;acs_im1(64:end,:,:)];
figure, imshow(rssq(acs_im),[]), title('CONC ACS Image')

acs_f = fft2(acs_im);
figure, imshow(rssq(acs_f),[]), title('CONC Fourier Domain')

%% DO the same for the MB Data
phases = exp(sqrt(-1)*pi*(0:n-1));
phas = repmat(phases,[m 1 no_c ims]);
nkspace = kspace;% .* phas;
phases2 = (exp(sqrt(-1)*(pi+(2*pi/m))*(0:m-1))).';
%phases2 = (exp(sqrt(-1)*(pi)*(0:m-1))).'; 
phas2 = repmat(phases2,[1 n no_c ims]);
mnkspace = nkspace .* phas2;

new_kspace = mnkspace;
% new_kspace = kspace;
figure, imshow(rssq(ifft2(new_kspace(:,:,:,1))),[]), title('MB Image for 1 out of 15')
figure, imshow(rssq(ifft2(new_kspace(:,:,:,10))),[]), title('MB Image for 10 out of 15')

save('inital_s_16','acs_f','new_kspace','slice_R','PE_R')
save('inital_s_16_ind','ACS1_final','ACS2_final','ACS3_final')

% coils = zeros(126,64,no_c,3);
% coils(:,:,:,1) = ACS1_final;
% coils(:,:,:,2) = ACS2_final;
% coils(:,:,:,3) = ACS3_final;
% save('coils_11','coils')