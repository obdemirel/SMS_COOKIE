clear all

save_num = 2;



save_name = ['cine_r1_reg_paper_' num2str(save_num)];
save_num = save_num+1;

load pre_data_golden.mat
load pre_data_ind_golden.mat

for asd3 = 1:4
dum1 = reshape(squeeze(data_ak_ind(asd3,:,:)),kernel_s_ind,kernel_s_ind,32,32);
        dum1 = permute(dum1,[2 1 3 4]);
    for asd1 = 1:32
        for asd2 = 1:32
        dum1(:,:,asd1,asd2) = ((dum1(:,:,asd1,asd2)).');
        end
    end
    data_ak_ind(asd3,:,:) = reshape(dum1,[kernel_s_ind*kernel_s_ind*32 32]);
end  

load('C:\Users\obura\Desktop\dMRI\spsg\sg_kspace_set1_recon.mat')
%     hf_kspace = permute(hf_kspace,[4 1 2 3 5]);
hf_kspace = permute(hf_kspace,[2 3 4 1 5]);
load('dmri_sense_maps_48_48acs_thr0_05.mat')
Coil_sensitivites = sense_maps;


DTDKMB = zeros(4*size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),size(hf_kspace,5),'single');
lambda_orig_s = zeros(4*size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),size(hf_kspace,5),'single');
ATb = zeros(4*size(hf_kspace,1)*size(hf_kspace,2)*size(hf_kspace,3),size(hf_kspace,5),'single');
x = zeros(4*size(hf_kspace,1)*size(hf_kspace,2)*size(hf_kspace,3),size(hf_kspace,5),'single');
sense_images =zeros(size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,5),4,'single');
pp_all = zeros(size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),4,size(hf_kspace,5),'single');

% data_kspace = data_kspace./max(data_kspace(:));
% kspace_grappa_recons = kspace_grappa_recons./max(kspace_grappa_recons(:));
cg_iter = 5;
outer_loop = 8; %10+1;

slices = 4;
ksb = size(hf_kspace,1);
num_images = size(hf_kspace,5);


lambda = 1e-1; %% should be 0

sigmas = [1 1 1 1];
rho = 1e-1;
llr_th_sca = 0.2; % 4.4 and 4.7
p1 = 8;
soft_sign = 1; %% zero means hard

before_c = 1;
after_c = 1;
%%% DO NOT FORGET FOR THE COOKIE%%% DO NOT FORGET FOR THE COOKIE
%%% DO NOT FORGET FOR THE COOKIE%%% DO NOT FORGET FOR THE COOKIE
%%% preperation of x z and K%%% DO NOT FORGET FOR THE COOKIE
for ss=1:num_images%%% DO NOT FORGET FOR THE COOKIE
    %%% fix the phases (incoming slice-MB is phase shifted)%%%%%%%%%%%%%%%%
    orig_slice = hf_kspace(:,:,:,:,ss);  %%% DO NOT FORGET FOR THE COOKIE
    orig_s(1:ksb,:,:)  = orig_slice(:,:,:,1);  %% original k-space
    orig_s(ksb*1 + 1:2*ksb,:,:)  = orig_slice(:,:,:,2); %% from slice-grappa
    orig_s(ksb*2 + 1:3*ksb,:,:)  = orig_slice(:,:,:,3);
    orig_s(ksb*3 + 1:4*ksb,:,:)  = orig_slice(:,:,:,4);
    orig_s = orig_s.*1.66;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kmb = (squeeze(data_kspace(ss,1:slices:end,:,:)));
    input = [kmb;kmb;kmb;kmb];
    
    [m,n,no_c] = size(input);
    
    %%% sampling points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    non_acq_p = kmb==0;
    acq_p = ones(ksb,n,no_c,'single')-non_acq_p;
    non_acq_p = logical(non_acq_p);
    loc_mask = logical(acq_p);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% some usefull functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);
    kspace_to_im = @(x) ifft2(x);
    im_to_kspace = @(x) fft2(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = @(x) encoder(x,m,n,no_c,slices,Coil_sensitivites,ksb);
    ET = @(x) encoder_t(x,m,n,no_c,slices,Coil_sensitivites,ksb);
    
    
    %% the whole ATA matrix is inside
    mat_op = @(x) matrix_maker(x,m,n,no_c,data_ak_ind,kernel_s_ind,ksb,slices,lambda,sigmas,rho,loc_mask,Coil_sensitivites);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
    DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% DT(D(KMB)) operation equals to select
    select = reshape(DT(D(kmb)),[ksb n no_c]);
    DTDKMB(:,:,:,ss) = [select;select;select;select];
    
    lambda_orig_s(:,:,:,ss) = lambda.*orig_s;
    
    
    %orig_s = zeros(m,n,no_c,'single');
    z(:,ss) = reshape(E(input(:)),m*n,1);
    %K(:,ss) = rho*z(:,ss);
    K = zeros(size(z));
    
    
    ATA =  @(x) mat_op(x);
    
    lls = lambda_orig_s(:,:,:,ss);
    ddts = DTDKMB(:,:,:,ss);
    
    
    %ATb = @(a,b,c,d,e) atb_func(a,b,c,d,e);
    %ATb(:,ss) = lls(:) + ddts(:) + 0.5*rho*(z(:,ss) - K(:,ss)/rho);
    if(lambda==0)
        %x(:,ss) = 1.*ddts(:);  %% should be 0;
        x(:,ss) = ddts(:);
    else
        x(:,ss) = ddts(:);
    end
end



for main_loop = 1:outer_loop
    
    for ss=1:num_images
        tic
        
        disp(['outer loop at: ' num2str(main_loop) ', T1 image of: ' num2str(ss)])
        tic
        lls = lambda_orig_s(:,:,:,ss);
        ddts = DTDKMB(:,:,:,ss);
        ATb = @(zzz,kkk) atb_func(ddts,lls,zzz,kkk,rho,Coil_sensitivites,ksb);
        if(main_loop==1)
%             [x(:,ss),error] = conjgrad(5,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss));
            [x(:,ss)] = conjgrad_clean(cg_iter,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss));
        else
            %[x(:,ss),error] = conjgrad(cg_iter,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss));
            [x(:,ss)] = conjgrad_clean(cg_iter,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss));
        end
        
        
        %%% z iteration
        z_plus_k(:,:,ss) = E(x(:,ss)) + reshape(K(:,ss)/rho,m,n);
        %z_image = reshape(z(:,ss),m,n,no_c);
        
        toc
    end
    
    tic
    %before_x = sense_images;
    
    %%%% before images
    for ss =1:num_images
        imbe1 = reshape(x(:,ss),[m n no_c]);
        image_before1(ss,:,:) = rssq(kspace_to_im((imbe1(1:ksb,:,:))));
        image_before2(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*1 + 1:2*ksb,:,:))));
        image_before3(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*2 + 1:3*ksb,:,:))));
        image_before4(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*3 + 1:4*ksb,:,:))));
    end
    image_before(before_c,1,:,:,:) = image_before1;
    image_before(before_c,2,:,:,:) = image_before2;
    image_before(before_c,3,:,:,:) = image_before3;
    image_before(before_c,4,:,:,:) = image_before4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    before_c = before_c+1;
    
    
    
    %% slice by slice accross all 15 images low rank denoising
    slice_n = slices;
    for slice_no =1:slice_n
        sense1 = z_plus_k((slice_no-1)*ksb + 1:slice_no*ksb,:,:);
        llr_th = llr_th_sca* max(abs(sense1(:)));
        denoised_sense1(:,:,:,slice_no) = llr_threshold(sense1, p1, p1, llr_th, soft_sign);
    end
    
    %%% after images
                        for ss =1:num_images
                            image_after1(ss,:,:) = denoised_sense1(:,:,ss,1);
                            image_after2(ss,:,:) = denoised_sense1(:,:,ss,2);
                            image_after3(ss,:,:) = denoised_sense1(:,:,ss,3);
                            image_after4(ss,:,:) = denoised_sense1(:,:,ss,4);
                        end
                        image_after(after_c,1,:,:,:) = image_after1;
                        image_after(after_c,2,:,:,:) = image_after2;
                        image_after(after_c,3,:,:,:) = image_after3;
                        image_after(after_c,4,:,:,:) = image_after4;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        after_c = after_c+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    
     for ss = 1:num_images
        %%%time to go back to k-space
        dd1 = [squeeze(denoised_sense1(:,:,ss,1));squeeze(denoised_sense1(:,:,ss,2));squeeze(denoised_sense1(:,:,ss,3));squeeze(denoised_sense1(:,:,ss,4))];
        z(:,ss) = dd1(:);
        xx = E(x(:,ss));
        K(:,ss) = K(:,ss) + rho*(xx(:)-z(:,ss));
        
    end
    %% error term
    %     resi1(main_loop) = norm(x(:,1)-z(:,1))/norm(x(:,1));
    %     disp(['Residual for 1: ' num2str(resi1(main_loop))])
    %     resi5(main_loop) = norm(x(:,5)-z(:,5))/norm(x(:,5));
    %     disp(['Residual for 5: ' num2str(resi5(main_loop))])
    %     resi10(main_loop) = norm(x(:,10)-z(:,10))/norm(x(:,10));
    %     disp(['Residual for 10: ' num2str(resi10(main_loop))])
    
    disp('Denoising...')
    toc
    %name = ['it_' num2str(main_loop)];
    %save(name,'x','z');
    
    %% plot the final results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       figure, plot(log(error(2:end)))
    %if(mod(main_loop,5)==1)
                        asd1 = abs(squeeze(image_before(after_c-1,3,1,:,:)));
                        asd1 = asd1./max(asd1(:));
                        asd2 = abs(squeeze(image_after(before_c-1,3,1,:,:)));
                        asd2 = asd2./max(asd2(:));
                        figure, imshow(flipud(cat(2,asd1,asd2)),[]), drawnow()
                        figure, imshow(flipud(abs(asd1-asd2)),[]), colorbar
                        
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for ss=1:num_images
    final_kspaces(ss,:,:,:) = reshape(x(:,ss),[m n no_c]);
end

[ims,m,n,no_c] = size(final_kspaces);
kspace_final = permute(final_kspaces,[2 3 4 1]);


[m,n,no_c,ims] = size(kspace_final);
mod_kspace = zeros(m/slices,n,no_c,slices,ims,'single');
ksb = size(hf_kspace,1);

ind = 1;
for im=1:ims
    for slice_num = [1 2 3 4]
        dummy_im2 = ifft2(kspace_final((slice_num-1)*ksb +1:slice_num*ksb,:,:,im));
        mod_kspace(:,:,:,slice_num,im) = fft2(dummy_im2);
    end
end

%     kspace_grappa_recons = mod_kspace;
%     kspace_recon_center = zeros(size(mod_kspace));
%     kspace_recon_center(:,91-16:91+16-1,:,:,:) = mod_kspace(:,91-16:91+16-1,:,:,:);
%
%
%
%     [img_spirit_sense1, img_spirit_phase_sens] = generate_images_MB_original(kspace_grappa_recons, kspace_recon_center);
%
%
save(save_name,'mod_kspace','lambda','sigmas','rho','llr_th_sca','-v7.3')
save('mid_images_llr','image_before','image_after','-v7.3')

%save('l4_s3_r2125e4_th5e2_8_soft_images','image_after','image_before','lambda','sigmas','rho','p1','llr_th_sca','soft_sign')

