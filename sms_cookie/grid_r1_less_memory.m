clear all

save_num = 1;

lambdas = [0];
sigmas_s  = [1];

for lambda_loop  = 1:size(lambdas,2)
    for simga_loop = 1:size(sigmas_s,2)
        
        save_name = ['dmri_wo_reg_deneme_' num2str(save_num)];
        save_num = save_num+1;

% load('C:\Users\obura\Desktop\dMRI\rock\initialize\pre_data_golden.mat')
load pre_data
load('pre_data_ind_golden.mat')

for asd3 = 1:4
dum1 = reshape(squeeze(data_ak_ind(asd3,:,:)),5,5,32,32);
        dum1 = permute(dum1,[2 1 3 4]);
    for asd1 = 1:32
        for asd2 = 1:32
        dum1(:,:,asd1,asd2) = ((dum1(:,:,asd1,asd2)).');
        end
    end
    data_ak_ind(asd3,:,:) = reshape(dum1,[5*5*32 32]);
end
    

slices = 4;
ksb = 234;

lambda = lambdas(lambda_loop);
sigmas = [sigmas_s(simga_loop) sigmas_s(simga_loop) sigmas_s(simga_loop) sigmas_s(simga_loop)];

[ims,m,n,no_c] = size(data_kspace);
mod_kspace = zeros(m/4,n,no_c,4,ims,'single');

for ss=1:5
    disp([num2str(ss)])
    tic
    
    %%% fix the phases (incoming slice-MB is phase shifted)%%%%%%%%%%%%%%%%
    
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

    input = [kmb;kmb;kmb;kmb];
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
    DTDKMB = [select;select;select;select];
    
    lambda_orig_s = lambda.*DTDKMB;
    
    ATA =  @(x) mat_op(x);
    ATb = lambda_orig_s + DTDKMB;
    
%     x = 0.*orig_s;
    x = zeros(size(DTDKMB),'single');
    
    %% cg iterations
    [x,error] = conjgrad(40,ATA, ATb(:), x(:));
    
    
    final_kspaces = reshape(x,[m n no_c]);
  
    mod_kspace(:,:,:,1,ss) = final_kspaces(1:234,:,:);
    mod_kspace(:,:,:,2,ss) = final_kspaces(235:2*234,:,:);
    mod_kspace(:,:,:,3,ss) = final_kspaces(469:3*234,:,:);
    mod_kspace(:,:,:,4,ss) = final_kspaces(703:4*234,:,:);
    
    %%%%% plot the final results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         a1 = rssq(ifft2(final_kspace(1:234,:,:))).';
%         a2 = rssq(ifft2(final_kspace(235:2*234,:,:))).';
%         a3 = rssq(ifft2(final_kspace(469:3*234,:,:))).';
%         a4 = rssq(ifft2(final_kspace(703:4*234,:,:))).';
%         a = cat(2,a1,a2,a3,a4);
%         figure, imshow(a,[]);
%         figure, plot(log(error(2:end)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    toc
end

save(save_name,'mod_kspace','lambdas','sigmas','-v7.3')
        
    end
end