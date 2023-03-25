clear all

save_num = 2;



save_name = ['cine_sms_reg__' num2str(save_num)];
save_num = save_num+1;


load pre_data_golden.mat
%data_kspace = squeeze(data_kspace(1:5,:,:,:));
data_kspace = (data_kspace(:,:,:,:));
load pre_data_ind_golden.mat
load sense_grappa_res


load sense_maps
Coil_sensitivites = sense_maps;

% 
% DTDKMB = zeros(5*size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),size(data_kspace,1),'single');
% lambda_orig_s = zeros(5*size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),size(data_kspace,1),'single');
% ATb = zeros(5*size(hf_kspace,1)*size(hf_kspace,2)*size(hf_kspace,3),size(data_kspace,1),'single');
% x = zeros(5*size(hf_kspace,1)*size(hf_kspace,2)*size(hf_kspace,3),size(data_kspace,1),'single');
% sense_images =zeros(size(hf_kspace,1),size(hf_kspace,2),size(data_kspace,1),3,'single');
% pp_all = zeros(size(hf_kspace,1),size(hf_kspace,2),size(hf_kspace,3),3,size(data_kspace,1),'single');

% data_kspace = data_kspace./max(data_kspace(:));
% kspace_grappa_recons = kspace_grappa_recons./max(kspace_grappa_recons(:));
cg_iter = 5;
outer_loop = 10; %10+1;

slices = 3;
ksb = size(data_kspace,2)/slices;
num_images = size(data_kspace,1);


lambda = 0; %% should be 0

sigmas = [1 1 1];% 1 1];
rho = 1e-1;
llr_th_sca = 0.075; % 4.4 and 4.7
p1 = 8;
soft_sign = 1; %% zero means hard

% rho = 1e-2;
% threshold_factor = 0.02; 

before_c = 1;
after_c = 1;
%%% DO NOT FORGET FOR THE COOKIE%%% DO NOT FORGET FOR THE COOKIE
%%% DO NOT FORGET FOR THE COOKIE%%% DO NOT FORGET FOR THE COOKIE
%%% preperation of x z and K%%% DO NOT FORGET FOR THE COOKIE
for ss=1:num_images%%% DO NOT FORGET FOR THE COOKIE
    %%% fix the phases (incoming slice-MB is phase shifted)%%%%%%%%%%%%%%%%
%     orig_slice = hf_kspace(:,:,:,:,ss);  %%% DO NOT FORGET FOR THE COOKIE
%     orig_s(1:ksb,:,:)  = orig_slice(:,:,:,1);  %% original k-space
%     orig_s(ksb*1 + 1:2*ksb,:,:)  = orig_slice(:,:,:,2); %% from slice-grappa
%     orig_s(ksb*2 + 1:3*ksb,:,:)  = orig_slice(:,:,:,3);
%     orig_s = orig_s.*1.66;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kmb = (squeeze(data_kspace(ss,1:slices:end,:,:)));
    input = [kmb;kmb;kmb];%;kmb;kmb];
    
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
    DTDKMB(:,:,:,ss) = [select;select;select];%;select;select];
    
    lambda_orig_s(:,:,:,ss) = lambda.*input;
    
    
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
        x(:,ss) = 0.*ddts(:);
        selec_final_kspaces = final_kspaces(:,:,:,ss);
		selec_final_kspaces = ifft2(selec_final_kspaces);
		a1 = fft2(selec_final_kspaces(1:348,:,:));
		a2 = fft2(selec_final_kspaces((1*348) + 1:2*348,:,:));
		a3 = fft2(selec_final_kspaces((2*348) + 1:3*348,:,:));
		%a4 = fft2(final_kspaces((3*348) + 1:4*348,:,:));
		%a5 = fft2(final_kspaces((4*348) + 1:5*348,:,:));
		
		selec_final_kspaces = [a1;a2;a3];%;a4;a5];
		
		diff_btw = max(abs(selec_final_kspaces(:)))./max(abs(kmb(:)));
		x(:,ss) = selec_final_kspaces(:)./diff_btw;
        z(:,ss) = reshape(E(x(:,ss)),m*n,1);
    else
        x(:,ss) = ddts(:);
    end
end


parpool(24)
for main_loop = 1:outer_loop
    
    parfor ss=1:num_images
        tic
        
        disp(['outer loop at: ' num2str(main_loop) ', T1 image of: ' num2str(ss)])
        tic
        lls = lambda_orig_s(:,:,:,ss);
        ddts = DTDKMB(:,:,:,ss);
        ATb = @(zzz,kkk) atb_func(ddts,lls,zzz,kkk,rho,Coil_sensitivites,ksb,slices);
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
        %image_before4(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*3 + 1:4*ksb,:,:))));
        %image_before5(ss,:,:) = rssq(kspace_to_im((imbe1(ksb*4 + 1:5*ksb,:,:))));
    end
    image_before(before_c,1,:,:,:) = image_before1;
    image_before(before_c,2,:,:,:) = image_before2;
    image_before(before_c,3,:,:,:) = image_before3;
    %image_before(before_c,4,:,:,:) = image_before4;
    %image_before(before_c,5,:,:,:) = image_before5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    before_c = before_c+1;
    
    
    
    slice_n = slices;
    for slice_no=1:slice_n
        %% Open for LOST
        %         x_coils(:,:,1) =  z_plus_k(0*ksb/slices + 1:1*ksb/slices,:,ss);
        %         x_coils(:,:,2) =  z_plus_k(1*ksb/slices + 1:2*ksb/slices,:,ss);
        %         x_coils(:,:,3) =  z_plus_k(2*ksb/slices + 1:3*ksb/slices,:,ss);
        %         x_coils(:,:,4) =  z_plus_k(3*ksb/slices + 1:4*ksb/slices,:,ss);
        %         x_coils(:,:,5) =  z_plus_k(4*ksb/slices + 1:5*ksb/slices,:,ss);
        %         [x_coils_alt] = LOST3(x_coils,threshold_factor);
        %         denoised_sense(:,:,:,ss) = x_coils_alt;
        sense1 = z_plus_k((slice_no-1)*ksb + 1:slice_no*ksb,:,:);
        
        if(slice_no==1)
            sense1 = circshift(sense1,[0 0 0]);
        elseif(slice_no==2)
            sense1 = circshift(sense1,[0 45 0]);
        elseif(slice_no==3)
            sense1 = circshift(sense1,[0 90 0]);
        %elseif(slice_no==4)
        %    sense1 = circshift(sense1,[0 135 0]);
        %elseif(slice_no==5)
        %    sense1 = circshift(sense1,[0 180 0]);
        end
        
        llr_th = llr_th_sca* max(abs(sense1(:)));
        denoised_sense1(:,:,:,slice_no) = llr_threshold(sense1, p1, p1, llr_th, soft_sign);
        
    end
    %denoised_sense1 = permute(denoised_sense,[1 2 4 3]); %% open for LOST
    denoised_sense1(:,:,:,1) = circshift(denoised_sense1(:,:,:,1),[0 0 0]);
    denoised_sense1(:,:,:,2) = circshift(denoised_sense1(:,:,:,2),[0 -45 0]);
    denoised_sense1(:,:,:,3) = circshift(denoised_sense1(:,:,:,3),[0 -90 0]);
    %denoised_sense1(:,:,:,4) = circshift(denoised_sense1(:,:,:,4),[0 -135 0]);
    %denoised_sense1(:,:,:,5) = circshift(denoised_sense1(:,:,:,5),[0 -180 0]);
    

    %%% after images
                        for ss =1:num_images
                            image_after1(ss,:,:) = denoised_sense1(:,:,ss,1);
                            image_after2(ss,:,:) = denoised_sense1(:,:,ss,2);
                            image_after3(ss,:,:) = denoised_sense1(:,:,ss,3);
                            %image_after4(ss,:,:) = denoised_sense1(:,:,ss,4);
                            %image_after5(ss,:,:) = denoised_sense1(:,:,ss,5);
                        end
                        image_after(after_c,1,:,:,:) = image_after1;
                        image_after(after_c,2,:,:,:) = image_after2;
                        image_after(after_c,3,:,:,:) = image_after3;
                        %image_after(after_c,4,:,:,:) = image_after4;
                        %image_after(after_c,5,:,:,:) = image_after5;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        after_c = after_c+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    
     for ss = 1:num_images
        %%%time to go back to k-space
        dd1 = [squeeze(denoised_sense1(:,:,ss,1));squeeze(denoised_sense1(:,:,ss,2));squeeze(denoised_sense1(:,:,ss,3))];%;squeeze(denoised_sense1(:,:,ss,4));squeeze(denoised_sense1(:,:,ss,5))];
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
    %                    asd1 = abs(squeeze(image_before(after_c-1,2,1,:,:)));
    %                    asd1 = asd1./max(asd1(:));
    %                    asd2 = abs(squeeze(image_after(before_c-1,2,1,:,:)));
    %                    asd2 = asd2./max(asd2(:));
    %                    figure, imshow(flipud(cat(2,asd1,asd2)),[]), drawnow()
    %                    figure, imshow(flipud(abs(asd1-asd2)),[]), colorbar
                        
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for ss=1:num_images
    final_kspaces2(ss,:,:,:) = reshape(x(:,ss),[m n no_c]);
end

[ims,m,n,no_c] = size(final_kspaces2);
kspace_final = permute(final_kspaces2,[2 3 4 1]);


[m,n,no_c,ims] = size(kspace_final);
mod_kspace = zeros(m/slices,n,no_c,slices,ims,'single');



ind = 1;
for im=1:ims
    for slice_num = [1 2 3];% 4 5]
        dummy_im2 = ifft2(kspace_final((slice_num-1)*ksb +1:slice_num*ksb,:,:,im));
        mod_kspace(:,:,:,slice_num,im) = fft2(dummy_im2);
    end
end

%save(save_name,'mod_kspace','lambda','rho','llr_th_sca','-v7.3')
% save(save_name,'mod_kspace','lambda','rho','threshold_factor','-v7.3')
image_after_save = squeeze(image_after(:,2,:,:,:));
image_before_save = squeeze(image_before(:,2,:,:,:));
save('mid_images_lost','image_after_save','image_before_save','rho','llr_th_sca','-v7.3')
% sense_maker

% sense_maps = Coil_sensitivites;
kspace_to_im = @(x) ifft2(x) * sqrt(size(x,1) * size(x,2));
sense1_images = zeros(348,178,3,23,'single');

for asd1 = 1:size(mod_kspace,5)
asd1
a1 = squeeze(mod_kspace(:,:,:,1,asd1));
a2 = squeeze(mod_kspace(:,:,:,2,asd1));
a3 = squeeze(mod_kspace(:,:,:,3,asd1));
%a4 = squeeze(mod_kspace(:,:,:,4,asd1));
%a5 = squeeze(mod_kspace(:,:,:,5,asd1));

a1 = squeeze((sum(conj(squeeze(sense_maps(:,:,:,1))) .* kspace_to_im(a1),3)));
a2 = squeeze((sum(conj(squeeze(sense_maps(:,:,:,2))) .* kspace_to_im(a2),3)));
a3 = squeeze((sum(conj(squeeze(sense_maps(:,:,:,3))) .* kspace_to_im(a3),3)));
%a4 = squeeze((sum(conj(squeeze(sense_maps(:,:,:,4))) .* kspace_to_im(a4),3)));
%a5 = squeeze((sum(conj(squeeze(sense_maps(:,:,:,5))) .* kspace_to_im(a5),3)));

sense1_images(:,:,1,asd1) = a1;
sense1_images(:,:,2,asd1) = a2;
sense1_images(:,:,3,asd1) = a3;
%sense1_images(:,:,4,asd1) = a4;
%sense1_images(:,:,5,asd1) = a5;
end


save('sms_llr_images','sense1_images','rho','llr_th_sca');

delete(gcp('nocreate'))