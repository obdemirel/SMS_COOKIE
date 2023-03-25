clear all

save_num = 1;

lambdas = [0];
sigmas_s  = [1];

for lambda_loop  = 1:size(lambdas,2)
    for simga_loop = 1:size(sigmas_s,2)
        
        save_name = ['dmri_wo_reg_deneme_5_' num2str(save_num)];
        save_num = save_num+1;
        
        
        
% load('C:\Users\obura\Desktop\dMRI\rock\initialize\pre_data_golden.mat')
load pre_data_golden
dummy = data_kspace;
dummy(:,:,1:2:end,:,:) = 0;
dummy(:,:,77:100,:,:) = data_kspace(:,:,77:100,:,:);
data_kspace = dummy;
load('pre_data_ind_golden.mat')
load sense_grappa_res


slices = 3;
ksb = 348;

lambda = lambdas(lambda_loop);
sigmas = [sigmas_s(simga_loop) sigmas_s(simga_loop) sigmas_s(simga_loop)];% sigmas_s(simga_loop) sigmas_s(simga_loop)];

parpool(12)
parfor ss=1:23
    disp([num2str(ss)])
    tic
    
    %%% fix the phases (incoming slice-MB is phase shifted)%%%%%%%%%%%%%%%%
    
%     orig_slice = hf_kspace(:,:,:,:,ss);
%     [m,n,no_c] = size(orig_slice(:,:,:,1));
%     phases2 = (exp(sqrt(-1)*(pi+(2*pi/m))*(0:m-1))).';
%     phas2 = 1;%repmat(phases2,[1 n no_c]);
%     orig_s(1:ksb,:,:)  = orig_slice(:,:,:,1) .* phas2;  %% original k-space
%     orig_s(ksb*1 + 1:2*ksb,:,:)  = orig_slice(:,:,:,2) .* phas2; %% from slice-grappa
%     orig_s(ksb*2 + 1:3*ksb,:,:)  = orig_slice(:,:,:,3) .* phas2;
%     orig_s(ksb*3 + 1:4*ksb,:,:)  = orig_slice(:,:,:,4) .* phas2;
%     
%     orig_s = orig_s.*1.66;
%     
%     [m,n,no_c] = size(orig_slice(:,:,:,1));
% 
%     orig_s(1:292,:,:)  = orig_slice(:,:,:,1);  %% original k-space
%     orig_s(293:2*292,:,:)  = orig_slice(:,:,:,2); %% from slice-grappa
%     orig_s(585:3*292,:,:)  = orig_slice(:,:,:,3);
%     orig_s(877:4*292,:,:)  = orig_slice(:,:,:,4);
%     orig_s(1169:5*292,:,:)  = orig_slice(:,:,:,5);
% %     orig_s = orig_s.*2;
    %orig_s = [squeeze(orig_slice(:,:,:,1));squeeze(orig_slice(:,:,:,2));squeeze(orig_slice(:,:,:,3))];
    %%%%orig_s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kmb = squeeze(data_kspace(ss,1:slices:end,:,:));
    [m,n,no_c] = size(kmb);
    %%% sampling points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    non_acq_p = kmb==0;
    acq_p = ones(ksb,n,no_c,'single')-non_acq_p;
    non_acq_p = logical(non_acq_p);
    loc_mask = logical(acq_p);
    %loc_mask = circshift(loc_mask, [0 -shifter(ss) 0 ]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    input = [kmb;kmb;kmb];%;kmb;kmb];
    [m,n,no_c] = size(input);
    

    %% some usefull functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);
    kspace_to_im = @(x) ifft2(x);
    im_to_kspace = @(x) fft2(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %GMB = @(x) gmb_to_ks(data_ak_set,x,kernel_s,kernel_s,m,n,no_c,0);
    
    %% the whole ATA matrix is inside
    mat_op = @(x) matrix_maker(x,m,n,no_c,data_ak_ind,kernel_s_ind,ksb,slices,lambda,sigmas,loc_mask);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
    DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);
    
    %%% DT(D(KMB)) operation equals to select
    select = reshape(DT(D(kmb)),[ksb n no_c]);
    DTDKMB = [select;select;select];%;select;select];
    
    lambda_orig_s = lambda.*DTDKMB;
    
    ATA =  @(x) mat_op(x);
    ATb = lambda_orig_s + DTDKMB;
    
%     x = 0.*orig_s;
    x = zeros(size(DTDKMB),'single');
    
    selec_final_kspaces = ifft2(final_kspaces(:,:,:,ss));
    a1 = fft2(selec_final_kspaces(1:ksb,:,:));
    a2 = fft2(selec_final_kspaces((1*ksb) + 1:2*ksb,:,:));
    a3 = fft2(selec_final_kspaces((2*ksb) + 1:3*ksb,:,:));
    %a4 = fft2(final_kspaces((3*348) + 1:4*348,:,:));
    %a5 = fft2(final_kspaces((4*348) + 1:5*348,:,:));
    
    selec_final_kspaces = [a1;a2;a3];%;a4;a5];
    diff_btw = max(abs(selec_final_kspaces(:)))./max(abs(kmb(:)));
    
    x = selec_final_kspaces(:)./diff_btw;
    %x = DTDKMB;
    
    
    %% cg iterations
    [x,error] = conjgrad(25,ATA, ATb(:), x(:));
    
    
    final_kspace2 = reshape(x,[m n no_c]);
  
    
    %%%%% plot the final results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    a1 = rssq(ifft2(final_kspace2(1:ksb,:,:))).';
    %    a2 = rssq(ifft2(final_kspace2((1*ksb) + 1:2*ksb,:,:))).';
    %    a3 = rssq(ifft2(final_kspace2((2*ksb) + 1:3*ksb,:,:))).';
    %    %a4 = rssq(ifft2(final_kspace2((3*ksb) + 1:4*ksb,:,:))).';
    %    %a5 = rssq(ifft2(final_kspace2((4*ksb) + 1:5*ksb,:,:))).';
    %    a = cat(2,a1,a2,a3);%,a4,a5);
    %    figure, imshow(a,[]);
    %    figure, plot(log(error(2:end)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    final_kspaces2(ss,:,:,:) = final_kspace2;
    
    toc
end

[ims,m,n,no_c] = size(final_kspaces2);
kspace_final = permute(final_kspaces2,[2 3 4 1]);


[m,n,no_c,ims] = size(kspace_final);
mod_kspace = zeros(m/slices,n,no_c,slices,ims,'single');


for im=1:ims
mod_kspace(:,:,:,1,im) = kspace_final(1:ksb,:,:,im);
mod_kspace(:,:,:,2,im) = kspace_final((1*ksb) + 1:2*ksb,:,:,im);
mod_kspace(:,:,:,3,im) = kspace_final((2*ksb) + 1:3*ksb,:,:,im);
%mod_kspace(:,:,:,4,im) = kspace_final((3*ksb) + 1:4*ksb,:,:,im);
%mod_kspace(:,:,:,5,im) = kspace_final((4*ksb) + 1:5*ksb,:,:,im);
end
  
save(save_name,'mod_kspace','lambdas','sigmas','-v7.3')
        
    end
end
delete(gcp('nocreate'))