function [recon_reg_images,recon_reg] = SMS_COOKIE_reg(llr,lambdas,sigmas,spsg_initial,data_kspace,sense_maps,data_ak_set,kernel_r,kernel_s,slice_R,outer_loop,cg_iter,gui_on,par_on)



DTDKMB = zeros(slice_R*size(spsg_initial,2),size(spsg_initial,3),size(spsg_initial,4),size(spsg_initial,5),'single');
lambda_orig_s = zeros(slice_R*size(spsg_initial,2),size(spsg_initial,3),size(spsg_initial,4),size(spsg_initial,5),'single');
ATb = zeros(slice_R*size(spsg_initial,2)*size(spsg_initial,3)*size(spsg_initial,4),size(spsg_initial,5),'single');
x = zeros(slice_R*size(spsg_initial,2)*size(spsg_initial,3)*size(spsg_initial,4),size(spsg_initial,5),'single');
sense_images =zeros(size(spsg_initial,2),size(spsg_initial,3),size(spsg_initial,5),4,'single');
pp_all = zeros(size(spsg_initial,2),size(spsg_initial,3),size(spsg_initial,4),4,size(spsg_initial,5),'single');



ksb = size(spsg_initial,2);
num_images = size(spsg_initial,5);

lambda = lambdas;
sigmas = ones(1,slice_R) * sigmas;

rho = 1e-1;%1e-1;
llr_th_sca = llr; %0.1; % 4.4 and 4.7
p1 = 8;
soft_sign = 1; %% zero means hard


before_c = 1;
after_c = 1;
%%% DO NOT FORGET FOR THE COOKIE%%% DO NOT FORGET FOR THE COOKIE
%%% DO NOT FORGET FOR THE COOKIE%%% DO NOT FORGET FOR THE COOKIE
%%% preperation of x z and K%%% DO NOT FORGET FOR THE COOKIE
for ss=1:num_images%%% DO NOT FORGET FOR THE COOKIE
    
    orig_slice = permute(spsg_initial(:,:,:,:,ss),[2 3 4 1]);
    [m,n,no_c] = size(orig_slice(:,:,:,1));
    for slis = 1:slice_R
        orig_s(ksb*(slis-1) + 1:slis*ksb,:,:)  = orig_slice(:,:,:,slis);
    end
    mult_fac = max(abs(data_kspace(:))) ./ max(abs(orig_s(:)));
    orig_s = orig_s.*1;%mult_fac;
    
    
    kmb = squeeze(data_kspace(1:slice_R:end,:,:,ss));
    [m,n,no_c] = size(kmb);
    
    
    input = [];
    for kmb_no = 1:slice_R
        input = [input;kmb];%;kmb;kmb];
    end
    
    [m,n,no_c] = size(input);
    
    %%% sampling points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    non_acq_p = kmb==0;
    acq_p = ones(ksb,n,no_c,'single')-non_acq_p;
    non_acq_p = logical(non_acq_p);
    loc_mask = logical(acq_p);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% some usefull functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);
    kspace_to_im = @(x) ifft2c(x);
    im_to_kspace = @(x) fft2c(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = @(x) encoder(x,m,n,no_c,slice_R,sense_maps,ksb);
    ET = @(x) encoder_t(x,m,n,no_c,slice_R,sense_maps,ksb);
    
    
    %% the whole ATA matrix is inside
    mat_op = @(x) matrix_maker(x,m,n,no_c,data_ak_set,kernel_r,kernel_s,ksb,slice_R,lambda,sigmas,rho,loc_mask,sense_maps);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
    DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% DT(D(KMB)) operation equals to select
    select = reshape(DT(D(kmb)),[ksb n no_c]);
    
    dum = [];
    for select_no = 1:slice_R
         dum = [dum;select];%;select;select];
    end
    DTDKMB(:,:,:,ss) = dum;
    
    lambda_orig_s(:,:,:,ss) = lambda.*orig_s;
    
    
    z(:,ss) = reshape(E(input(:)),m*n,1);
    K = zeros(size(z));
    
    
    ATA =  @(x) mat_op(x);
    
    lls = lambda_orig_s(:,:,:,ss);
    ddts = DTDKMB(:,:,:,ss);
    

    if(lambda==0)
        x(:,ss) = ddts(:)*0;
    else
        x(:,ss) = orig_s(:)*1;
        z(:,ss) = reshape(E(orig_s(:)*1),m*n,1);
        %x(:,ss) = ddts(:)*0;
        %z(:,ss) = reshape(E(ddts(:)*0),m*n,1);
    end
end



for main_loop = 1:outer_loop
    
    parfor ss=1:num_images
        tic
        
        disp(['outer loop at: ' num2str(main_loop) ', T1 image of: ' num2str(ss)])
        tic
        lls = lambda_orig_s(:,:,:,ss);
        ddts = DTDKMB(:,:,:,ss);
        ATb = @(zzz,kkk) atb_func(ddts,lls,zzz,kkk,rho,sense_maps,ksb,slice_R);
        if(main_loop==1)
            [x(:,ss)] = conjgrad(cg_iter,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss),gui_on);
        else
            [x(:,ss)] = conjgrad(cg_iter,ATA, ATb(z(:,ss),K(:,ss)), x(:,ss),gui_on);
        end
        
        
        %%% z iteration
        z_plus_k(:,:,ss) = E(x(:,ss)) + reshape(K(:,ss)/rho,m,n);
        %z_image = reshape(z(:,ss),m,n,no_c);
        
        toc
    end
    
    tic
    %before_x = sense_images;
    
%     %%%% before images
%     for ss =1:num_images
%         imbe1 = reshape(x(:,ss),[m n no_c]);
%         image_before1(ss,:,:) = rssq(kspace_to_im((imbe1(1:ksb,:,:))));
%         image_before2(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*1 + 1:2*ksb,:,:))));
%         image_before3(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*2 + 1:3*ksb,:,:))));
%     end
%     image_before(before_c,1,:,:,:) = image_before1;
%     image_before(before_c,2,:,:,:) = image_before2;
%     image_before(before_c,3,:,:,:) = image_before3;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     before_c = before_c+1;
%     
    
    
    %% slice by slice accross all 15 images low rank denoising
    slice_n = slice_R;
    shift_amounts = ([1:slice_R]-1)*96;
    %shift_amounts = [46,0,-46];
    shift_amounts = [46,0,-46,46];
    parfor slice_no =1:slice_n
        sense1 = z_plus_k((slice_no-1)*ksb + 1:slice_no*ksb,:,:);
        sense1 = circshift(sense1,[0 shift_amounts(slice_no) 0]);
        llr_th = llr_th_sca* max(abs(sense1(:)));
        denoised_sense1(:,:,:,slice_no) = llr_threshold(sense1, p1, p1, llr_th, soft_sign);
    end
    
    for slice_no=1:slice_n
    denoised_sense1(:,:,:,slice_no) = circshift(denoised_sense1(:,:,:,slice_no),[0 -shift_amounts(slice_no) 0]);
    end
    
    %%% after images
%     for ss =1:num_images
%         image_after1(ss,:,:) = denoised_sense1(:,:,ss,1);
%         image_after2(ss,:,:) = denoised_sense1(:,:,ss,2);
%         image_after3(ss,:,:) = denoised_sense1(:,:,ss,3);
%     end
%     image_after(after_c,1,:,:,:) = image_after1;
%     image_after(after_c,2,:,:,:) = image_after2;
%     image_after(after_c,3,:,:,:) = image_after3;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     after_c = after_c+1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    
    for ss = 1:num_images
        %%%time to go back to k-space
        dd1 = [];
        for slis = 1:slice_R
            dd1 = [dd1;squeeze(denoised_sense1(:,:,ss,slis))];
        end
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
%     asd1 = abs(squeeze(image_before(after_c-1,2,2,:,:)));
%     asd1 = asd1./max(asd1(:));
%     asd2 = abs(squeeze(image_after(before_c-1,2,2,:,:)));
%     asd2 = asd2./max(asd2(:));
%     figure, imshow((cat(2,asd1,asd2)),[]), drawnow()
%     figure, imshow((abs(asd1-asd2)),[]), colorbar
    
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%% Saving 
for ss=1:num_images
    recon_reg(:,:,:,ss) = reshape(x(:,ss),[m n no_c]);
end

recon_reg_images = permute(denoised_sense1,[1 2 4 3]);

end

