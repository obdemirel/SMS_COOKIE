clear all

save_num = 1;

selector = [4e-4 0.04];  %% selector = [1e-4 0.04;4e-4 0.04];

for selector_loop = 1:size(selector,1)
    
    save_name = ['sp_sg_mb3_r4_2mm_64acs_recon_reg_' num2str(save_num)];
    save_num = save_num+1;
    
    load pre_data
    load pre_data_ind
    load sp_sg_mb3_r4_2mm_64acs_recon
    hf_kspace = permute(hf_kspace,[2 3 4 1 5]);
    
    
    
    DTDKMB = zeros(3*size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),size(hf_kspace,5),'single');
    lambda_orig_s = zeros(3*size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),size(hf_kspace,5),'single');
    ATb = zeros(3*size(hf_kspace,1)*size(hf_kspace,2)*size(hf_kspace,3),size(hf_kspace,5),'single');
    x = zeros(3*size(hf_kspace,1)*size(hf_kspace,2)*size(hf_kspace,3),size(hf_kspace,5),'single');
    sense_images =zeros(size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,5),3,'single');
    pp_all = zeros(size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),3,size(hf_kspace,5),'single');
    
    % data_kspace = data_kspace./max(data_kspace(:));
    % kspace_grappa_recons = kspace_grappa_recons./max(kspace_grappa_recons(:));
    cg_iter = 5;
    outer_loop = 25; %10+1;
    
    slices = 3;
    ksb = size(hf_kspace,1);
    num_images = size(hf_kspace,5);
    

    lambda = 1e-4; %% should be 0

    sigmas = [1e-3 1e-3 1e-3];
    rho = selector(selector_loop,1); %5.5e-4;
    llr_th_sca = selector(selector_loop,2); % 4.4 and 4.7
    p1 = 8;
    soft_sign = 1; %% zero means hard
    
    before_c = 1;
    after_c = 1;
    
    
    %%% preperation of x z and K
    for ss=1:num_images
        %%% fix the phases (incoming slice-MB is phase shifted)%%%%%%%%%%%%%%%%
        orig_slice = hf_kspace(:,:,:,:,ss);
        [m,n,no_c] = size(orig_slice(:,:,:,1));
        phases2 = (exp(sqrt(-1)*(pi+(2*pi/m))*(0:m-1))).';
        phas2 = 1;%repmat(phases2,[1 n no_c]);
        orig_s(1:ksb,:,:)  = orig_slice(:,:,:,1) .* phas2;  %% original k-space
        orig_s(ksb*1 + 1:2*ksb,:,:)  = orig_slice(:,:,:,2) .* phas2; %% from slice-grappa
        orig_s(ksb*2 + 1:3*ksb,:,:)  = orig_slice(:,:,:,3) .* phas2;
        %     orig_s = orig_s.*(1e-4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        kmb = (squeeze(data_kspace(ss,1:slices:end,:,:)));
        input = [kmb;kmb;kmb];
        
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
        
        %GMB = @(x) gmb_to_ks(data_ak_set,x,kernel_s,kernel_s,m,n,no_c,0);
        pp_all(:,:,:,1,ss) = orig_s(1:ksb,:,:);
        pp_all(:,:,:,2,ss) = orig_s(ksb*1 + 1:2*ksb,:,:);
        pp_all(:,:,:,3,ss) = orig_s(ksb*2 + 1:3*ksb,:,:);
        pp_all_center = zeros(size(pp_all));
        pp_all_center(center_locs(1):end,center_locs(2):center_locs(3),:,:,:) = pp_all(center_locs(1):end,center_locs(2):center_locs(3),:,:,:);
        [Coil_sensitivites] = generate_coils_sens(pp_all, pp_all_center);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        E = @(x) encoder(x,m,n,no_c,3,Coil_sensitivites,ksb);
        ET = @(x) encoder_t(x,m,n,no_c,3,Coil_sensitivites,ksb);
        
        
        %% the whole ATA matrix is inside
        mat_op = @(x) matrix_maker(x,m,n,no_c,data_ak_ind,kernel_s_ind,ksb,slices,lambda,sigmas,rho,loc_mask,Coil_sensitivites);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
        DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%% DT(D(KMB)) operation equals to select
        select = reshape(DT(D(kmb)),[ksb n no_c]);
        DTDKMB(:,:,:,ss) = [select;select;select];
        
        lambda_orig_s(:,:,:,ss) = lambda.*orig_s;
        
        
        
        z(:,ss) = reshape(E(orig_s(:)),m*n,1);
        %K(:,ss) = rho*z(:,ss);
        K = zeros(size(z));
        
        
        ATA =  @(x) mat_op(x);
        
        lls = lambda_orig_s(:,:,:,ss);
        ddts = DTDKMB(:,:,:,ss);
        
        
        %ATb = @(a,b,c,d,e) atb_func(a,b,c,d,e);
        %ATb(:,ss) = lls(:) + ddts(:) + 0.5*rho*(z(:,ss) - K(:,ss)/rho);
        

        x(:,ss) = 1.*ddts(:);

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
                [x(:,ss),error] = conjgrad(5,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss));
            else
                [x(:,ss),error] = conjgrad(cg_iter,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss));
            end
            
            
            %%% z iteration
            z_plus_k(:,:,ss) = E(x(:,ss)) + reshape(K(:,ss)/rho,m,n);
            %z_image = reshape(z(:,ss),m,n,no_c);
            
            toc
        end
        
        tic
        %before_x = sense_images;
        
        if(mod(main_loop,5)==1)
        %%%% before images
        for ss =1:num_images
            imbe1 = reshape(x(:,ss),[m n no_c]);
            image_before1(ss,:,:) = rssq(kspace_to_im((imbe1(1:ksb,:,:))));
            image_before2(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*1 + 1:2*ksb,:,:))));
            image_before3(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*2 + 1:3*ksb,:,:))));
        end
        image_before(before_c,1,:,:,:) = image_before1;
        image_before(before_c,2,:,:,:) = image_before2;
        image_before(before_c,3,:,:,:) = image_before3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        before_c = before_c+1;
        %mid_saver(x,m,n,no_c,15,main_loop)
        end
        
        
        %% slice by slice accross all 15 images low rank denoising
        slice_n = 3;
        for slice_no =1:slice_n
            sense1 = z_plus_k((slice_no-1)*ksb + 1:slice_no*ksb,:,:);
            llr_th = llr_th_sca* max(abs(sense1(:)));
            denoised_sense1(:,:,:,slice_no) = llr_threshold(sense1, p1, p1, llr_th, soft_sign);
            %after_x(:,:,:,slice_no) =denoised_sense1;
            %new_z((slice_no-1)*ksb + 1:slice_no*ksb,:,:,:) = permute(repmat(denoised_sense1,[1 1 1 no_c]),[1 2 4 3]).*repmat(coil_sens(:,:,:,slice_no),[1 1 1 15]);
            %new_z((slice_no-1)*ksb + 1:slice_no*ksb,:,:,:) = permute(repmat(denoised_sense1,[1 1 1 no_c]),[1 2 4 3]).*repmat(Coil_sensitivites(:,:,:,slice_no),[1 1 1 15]);
            %z_kspace((slice_no-1)*ksb + 1:slice_no*ksb,:,:,:) = im_to_kspace(new_z((slice_no-1)*ksb + 1:slice_no*ksb,:,:,:));
        end
        
        %%% after images
        if(mod(main_loop,5)==1)
        for ss =1:num_images
            image_after1(ss,:,:) = denoised_sense1(:,:,ss,1);
            image_after2(ss,:,:) = denoised_sense1(:,:,ss,2);
            image_after3(ss,:,:) = denoised_sense1(:,:,ss,3);
        end
        image_after(after_c,1,:,:,:) = image_after1;
        image_after(after_c,2,:,:,:) = image_after2;
        image_after(after_c,3,:,:,:) = image_after3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        after_c = after_c+1;
        %after_saver(z_kspace,m,n,no_c,15,main_loop)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        
        
        for ss = 1:num_images
            %%%time to go back to k-space
            dd1 = [denoised_sense1(:,:,ss,1);denoised_sense1(:,:,ss,2);denoised_sense1(:,:,ss,3)];
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
        
%         %% plot the final results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  %       figure, plot(log(error(2:end)))
%                 %if(mod(main_loop,5)==1)
%                     asd1 = abs(squeeze(image_before(after_c-1,1,1,:,:)));
%                     asd1 = asd1./max(asd1(:));
%         
%                     asd2 = abs(squeeze(image_after(before_c-1,1,1,:,:)));
%                     asd2 = asd2./max(asd2(:));
%         
%                     figure, imshow(cat(2,asd1,asd2),[])
%                     figure, imshow(abs(asd1-asd2),[]), colorbar
%         %end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    for ss=1:num_images
        final_kspaces(ss,:,:,:) = reshape(x(:,ss),[m n no_c]);
    end
    
    [ims,m,n,no_c] = size(final_kspaces);
    kspace_final = permute(final_kspaces,[2 3 4 1]);
    
    
    [m,n,no_c,ims] = size(kspace_final);
    mod_kspace = zeros(m/3,n,no_c,3,ims,'single');
    ksb = size(hf_kspace,1);
    
    ind = 1;
    for im=1:ims
        for slice_num = [1 2 3]
            dummy_im2 = ifft2(kspace_final((slice_num-1)*ksb +1:slice_num*ksb,:,:,im));
            mod_kspace(:,:,:,slice_num,im) = fft2(dummy_im2);
        end
    end
    
    kspace_grappa_recons = mod_kspace;
    kspace_recon_center = zeros(size(mod_kspace));
    kspace_recon_center(center_locs(1):end,center_locs(2):center_locs(3),:,:,:) = mod_kspace(center_locs(1):end,center_locs(2):center_locs(3),:,:,:);
    
    
    
    [img_spirit_sense1, img_spirit_phase_sens] = generate_images_MB_original(kspace_grappa_recons, kspace_recon_center);
    
    
    save(save_name,'img_spirit_phase_sens','lambda','sigmas','rho','p1','llr_th_sca','soft_sign','image_after','image_before')
    %save('l4_s3_r2125e4_th5e2_8_soft_images','image_after','image_before','lambda','sigmas','rho','p1','llr_th_sca','soft_sign')
    
    
end