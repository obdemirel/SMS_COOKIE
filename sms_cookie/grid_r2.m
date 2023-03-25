clear all

save_num = 1;

lambdas = [1e-4 0];
sigmas_s  = [1e-3];

for lambda_loop  = 1:size(lambdas,2)
    for simga_loop = 1:size(sigmas_s,2)
        
        save_name = ['cine_r2_wo_reg_paper_' num2str(save_num)];
        save_num = save_num+1;
        
        
load pre_data_cine
load initial_cine_kspace
load kspace_R2
%load res
% kspace_R1 = zeros(292,152,34,5,1);
% kspace_R1(:,:,:,:,1) = res;

slices = 5;
ksb = 292;

lambda = lambdas(lambda_loop);
sigmas = [sigmas_s(simga_loop) sigmas_s(simga_loop) sigmas_s(simga_loop) sigmas_s(simga_loop) sigmas_s(simga_loop)];


for ss=1:23
    disp([num2str(ss)])
    tic
    
    %%% fix the phases (incoming slice-MB is phase shifted)%%%%%%%%%%%%%%%%
    
%     orig_slice(:,:,:,1) = circshift(squeeze(kspace_MB5_R1(:,:,:,1,ss)), [31 -15 0 0]);
%     orig_slice(:,:,:,2) = circshift(squeeze(kspace_MB5_R1(:,:,:,2,ss)), [31 -15 0 0]);
%     orig_slice(:,:,:,3) = circshift(squeeze(kspace_MB5_R1(:,:,:,3,ss)), [30 -14 0 0]);
%     orig_slice(:,:,:,4) = circshift(squeeze(kspace_MB5_R1(:,:,:,4,ss)), [30 -14 0 0]);
%     orig_slice(:,:,:,5) = circshift(squeeze(kspace_MB5_R1(:,:,:,5,ss)), [30 -14 0 0]);
    orig_slice(:,:,:,1) = (squeeze(kspace_R2(:,:,:,1,ss)));
    orig_slice(:,:,:,2) = (squeeze(kspace_R2(:,:,:,2,ss)));
    orig_slice(:,:,:,3) = (squeeze(kspace_R2(:,:,:,3,ss)));
    orig_slice(:,:,:,4) = (squeeze(kspace_R2(:,:,:,4,ss)));
    orig_slice(:,:,:,5) = (squeeze(kspace_R2(:,:,:,5,ss)));
    
    [m,n,no_c] = size(orig_slice(:,:,:,1));

    orig_s(1:292,:,:)  = orig_slice(:,:,:,1);  %% original k-space
    orig_s(293:2*292,:,:)  = orig_slice(:,:,:,2); %% from slice-grappa
    orig_s(585:3*292,:,:)  = orig_slice(:,:,:,3);
    orig_s(877:4*292,:,:)  = orig_slice(:,:,:,4);
    orig_s(1169:5*292,:,:)  = orig_slice(:,:,:,5);
%     orig_s = orig_s.*2;
    %orig_s = [squeeze(orig_slice(:,:,:,1));squeeze(orig_slice(:,:,:,2));squeeze(orig_slice(:,:,:,3))];
    %%%%orig_s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kmb_full = squeeze(new_kspace(:,:,:,ss));
    kmb_r2 = zeros(size(kmb_full));
    kmb_r2(:,1:2:end,:) = kmb_full(:,1:2:end,:);
    kmb_r2(:,76-12:76+12-1,:) = kmb_full(:,76-12:76+12-1,:);
    kmb = circshift(kmb_r2, [0 0 0]);
    %%% sampling points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    non_acq_p = kmb==0;
    acq_p = ones(ksb,n,no_c,'single')-non_acq_p;
    non_acq_p = logical(non_acq_p);
    loc_mask = logical(acq_p);
    %loc_mask = circshift(loc_mask, [0 -shifter(ss) 0 ]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    input = [kmb;kmb;kmb;kmb;kmb];
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
    DTDKMB = [select;select;select;select;select];
    
    lambda_orig_s = lambda.*orig_s;
    
    ATA =  @(x) mat_op(x);
    ATb = lambda_orig_s + DTDKMB;
    
    x = 0.*orig_s;
    
    %% cg iterations
    [x,error] = conjgrad(25,ATA, ATb(:), x(:));
    
    
    final_kspace = reshape(x,[m n no_c]);
    
    %%%%% plot the final results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     a1 = rssq(ifft2(final_kspace(1:160,:,:)));
    %     a2 = rssq(ifft2(final_kspace(161:320,:,:)));
    %     a3 = rssq(ifft2(final_kspace(321:480,:,:)));
    %     a = cat(1,a1,a2,a3);
    %     figure, imshow(a,[]);
    %     figure, plot(log(error(2:end)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    final_kspaces(ss,:,:,:) = final_kspace;
    
    toc
end

[ims,m,n,no_c] = size(final_kspaces);
kspace_final = permute(final_kspaces,[2 3 4 1]);


[m,n,no_c,ims] = size(kspace_final);
mod_kspace = zeros(m/5,n,no_c,5,ims,'single');
ksb = 292;

for im=1:ims
mod_kspace(:,:,:,1,im) = kspace_final(1:292,:,:,im);
mod_kspace(:,:,:,2,im) = kspace_final(293:2*292,:,:,im);
mod_kspace(:,:,:,3,im) = kspace_final(585:3*292,:,:,im);
mod_kspace(:,:,:,4,im) = kspace_final(877:4*292,:,:,im);
mod_kspace(:,:,:,5,im) = kspace_final(1169:5*292,:,:,im);
end



kspace_grappa_recons = mod_kspace;
kspace_recon_center = zeros(size(kspace_grappa_recons));
kspace_recon_center(4:289,72-12:72+12-1,:,:,:) = mod_kspace(4:289,72-12:72+12-1,:,:,:);
[img_spirit_sense1, img_spirit_phase_sens] = generate_images_MB_original(kspace_grappa_recons, kspace_recon_center);
  
                
save(save_name,'img_spirit_sense1','lambda','sigmas')
        
    end
end