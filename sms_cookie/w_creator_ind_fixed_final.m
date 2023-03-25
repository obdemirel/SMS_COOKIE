%tic
clear all
%close all


load('C:\Users\obura\Desktop\dMRI\spsg\initial\inital_ind.mat')

for abc = 1:4

      if abc ==1
        conc_acs_fregion = ACS1_final;
    elseif abc == 2
        conc_acs_fregion = ACS2_final;
    elseif abc == 3
        conc_acs_fregion = ACS3_final;
        elseif abc == 4
        conc_acs_fregion = ACS4_final;
    end
    
    
    reg = 1e-2;
    
    %%% KERNEL %%%%%%%%%
    % XXX        XXX
    % XXX    ->  X1X
    % XXX        XXX
    %%%%%%%%%%%%%%%%%%%%
    
    kernel_row = 5;  %3
    kernel_col = 5;  %3
    
    kernel = calibSPIRiT(conc_acs_fregion, [kernel_row kernel_col], size(conc_acs_fregion,3), reg);
    ak = reshape(kernel,[size(kernel,1)*size(kernel,2)*size(kernel,3) size(kernel,4)]);
    
    
    
    data_ak_ind(abc,:,:) = ak;
    kernel_s_ind = kernel_row;
    
end


save('pre_data_ind_fixed','data_ak_ind','kernel_s_ind')

function kernel = calibSPIRiT(kCalib, kSize, nCoils, CalibTyk)

% kernel = calibSPIRiT(kCalib, kSize, nCoils, CalibTyk)
%
% Function calibrates a SPIRiT kernel from a calibration area in k-space
%
%
% (c) Michael Lustig 2013
%

[AtA] = dat2AtA(kCalib,kSize);
for n=1:nCoils
    tic
    kernel(:,:,:,n) = calibrate(AtA,kSize,nCoils,n,CalibTyk);
    disp(['Coil #: ',num2str(n),' is done out of ',num2str(nCoils)])
    toc
end
end

function [AtA,A,kernel] = dat2AtA(data, kSize)

% [AtA,A,kernel] = dat2AtA(data, kSize)
%
% Function computes the calibration matrix from calibration data.
%
% (c) Michael Lustig 2013



[sx,sy,nc] = size(data);

tmp = im2row(data,kSize); [tsx,tsy,tsz] = size(tmp);
A = reshape(tmp,tsx,tsy*tsz);

AtA = A'*A;

kernel = AtA;
kernel = reshape(kernel,kSize(1),kSize(2),nc,size(kernel,2));
end

function [kernel,rawkernel] = calibrate(AtA, kSize, nCoil, coil, lambda, sampling)

if nargin < 6
    sampling = ones([kSize,nCoil]);
end


dummyK = zeros(kSize(1),kSize(2),nCoil); dummyK((end+1)/2,(end+1)/2,coil) = 1;
idxY = find(dummyK);
sampling(idxY) = 0;
idxA = find(sampling);

Aty = AtA(:,idxY); Aty = Aty(idxA);
AtA = AtA(idxA,:); AtA =  AtA(:,idxA);

kernel = sampling*0;

lambda = norm(AtA,'fro')/size(AtA,1)*lambda;

rawkernel = inv(AtA + eye(size(AtA))*lambda)*Aty;
kernel(idxA) = rawkernel;
end

function res = im2row(im, winSize)
%res = im2row(im, winSize)
[sx,sy,sz] = size(im);

res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);
count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
            (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);
    end
end
end
