function [sense1_images,final_kspaces2] = SMS_COOKIE(lambdas,sigmas,spsg_initial,data_kspace,sense_maps,data_ak_set,kernel_r,kernel_s,slice_R,cg_iter,gui_on,par_on)




ksb = size(data_kspace,1)/slice_R;

lambda = lambdas;
sigmas = ones(1,slice_R) * sigmas;

if(par_on==1)
    parfor ss=1:size(data_kspace,4)
        tic
        
        orig_slice = permute(spsg_initial(:,:,:,:,ss),[2 3 4 1]);
        [m,n,no_c] = size(orig_slice(:,:,:,1));
        for slis = 1:slice_R
            orig_s(ksb*(slis-1) + 1:slis*ksb,:,:)  = orig_slice(:,:,:,slis);
        end
        
        if(lambda~=0)
        mult_fac = max(abs(data_kspace(:))) ./ max(abs(orig_s(:)));
        orig_s = orig_s.*1;
        end
        
        
        kmb = squeeze(data_kspace(1:slice_R:end,:,:,ss));
        [m,n,no_c] = size(kmb);
        
        %%% sampling points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        non_acq_p = kmb==0;
        acq_p = ones(ksb,n,no_c,'single')-non_acq_p;
        non_acq_p = logical(non_acq_p);
        loc_mask = logical(acq_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        input = [];
        for kmb_no = 1:slice_R
            input = [input;kmb];%;kmb;kmb];
        end
        [m,n,no_c] = size(input);
        
        
        
        
        %% the whole ATA matrix is inside
        mat_op = @(x) matrix_maker(x,m,n,no_c,data_ak_set,kernel_r,kernel_s,ksb,slice_R,lambda,sigmas,loc_mask);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
        DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);
        
        %%% DT(D(KMB)) operation equals to select
        select = reshape(DT(D(kmb)),[ksb n no_c]);
        DTDKMB = [];
        for select_no = 1:slice_R
            DTDKMB = [DTDKMB;select];%;select;select];
        end
        
        lambda_orig_s = lambda.*orig_s;
        %lambda_orig_s = 1.*orig_s;
        
        ATA =  @(x) mat_op(x);
        ATb = lambda_orig_s + DTDKMB;
        %ATb = lambda_orig_s + lambda.*DTDKMB;
        
        %x = zeros(size(DTDKMB),'single');
        if(lambda==0)
            x = DTDKMB*1;
        else
            x = orig_s*1;
        end
        
        
        
        %% cg iterations
        [x,error] = conjgrad(cg_iter,ATA, ATb(:), x(:),gui_on); %% conjugate gradient operations
        
        
        final_kspaces2(ss,:,:,:) = reshape(x,[m n no_c]);
        
        toc
    end
else
    for ss=1:size(data_kspace,4)
        tic
        
        orig_slice = permute(spsg_initial(:,:,:,:,ss),[2 3 4 1]);
        [m,n,no_c] = size(orig_slice(:,:,:,1));
        for slis = 1:slice_R
            orig_s(ksb*(slis-1) + 1:slis*ksb,:,:)  = orig_slice(:,:,:,slis);
        end
        
        if(lambda~=0)
        mult_fac = max(abs(data_kspace(:))) ./ max(abs(orig_s(:)));
        orig_s = orig_s.*1;%mult_fac;
        end
        
        kmb = squeeze(data_kspace(1:slice_R:end,:,:,ss));
        [m,n,no_c] = size(kmb);
        
        %%% sampling points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        non_acq_p = kmb==0;
        acq_p = ones(ksb,n,no_c,'single')-non_acq_p;
        non_acq_p = logical(non_acq_p);
        loc_mask = logical(acq_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        input = [];
        for kmb_no = 1:slice_R
            input = [input;kmb];%;kmb;kmb];
        end
        [m,n,no_c] = size(input);
        
        
        
        
        %% the whole ATA matrix is inside
        mat_op = @(x) matrix_maker(x,m,n,no_c,data_ak_set,kernel_r,kernel_s,ksb,slice_R,lambda,sigmas,loc_mask);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
        DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);
        
        %%% DT(D(KMB)) operation equals to select
        select = reshape(DT(D(kmb)),[ksb n no_c]);
        DTDKMB = [];
        for select_no = 1:slice_R
            DTDKMB = [DTDKMB;select];%;select;select];
        end
        
        lambda_orig_s = lambda.*orig_s;
        
        ATA =  @(x) mat_op(x);
        ATb = lambda_orig_s + DTDKMB;
        
        %x = zeros(size(DTDKMB),'single');
        if(lambda==0)
            x = DTDKMB*1;
        else
            x = orig_s;
        end
        
        
        
        %% cg iterations
        [x,error] = conjgrad(cg_iter,ATA, ATb(:), x(:),gui_on); %% conjugate gradient operations
        
        
        final_kspaces2(ss,:,:,:) = reshape(x,[m n no_c]);
        
        toc
        
        
    end
end


[ims,m,n,no_c] = size(final_kspaces2);
kspace_final = permute(final_kspaces2,[2 3 4 1]);

[m,n,no_c,ims] = size(kspace_final);
mod_kspace = zeros(m/slice_R,n,no_c,slice_R,ims,'single');


for im=1:ims
    for slis = 1:slice_R
        mod_kspace(:,:,:,slis,im) = kspace_final(((slis-1)*ksb) + 1:slis*ksb,:,:,im);
    end
end

ims = size(data_kspace,4);
[sense1_images] = sense1_maker(permute(mod_kspace,[4 1 2 3 5]),sense_maps,ims);


end

