function outp = matrix_maker(x,m,n,no_c,data_ak_ind,kernel_r_ind,kernel_s_ind,ksb,slices,lambda,sigmas,loc_mask)

x = reshape(x,m,n,no_c);
D = @(x) selection_operator(x,loc_mask,ksb,n,no_c); %locations specify
DT = @(x) adjoint_selection_operator(x,loc_mask,ksb,n,no_c);

for slis = 1:slices
    adder_all(:,:,:,slis) = DT(D(x(((slis-1)*ksb) + 1:ksb*slis,:,:)));
end

adder  =0;
for slis = 1:slices
adder = adder + reshape(squeeze(adder_all(:,:,:,slis)),[ksb n no_c]);
end

for abc = 1:slices
    
    ak_ind = squeeze(data_ak_ind(:,:,abc));
    sigma = sigmas(abc);
    
    x1 = x((abc-1)*ksb +1:abc*ksb,:,:);
    x1 = x1(:);
    
    
    G_ind = @(y) conv_op_ind(ak_ind,y,kernel_r_ind,kernel_s_ind,ksb,n,no_c);
    GT_ind = @(y) conv_op_ind_t(ak_ind,y,kernel_r_ind,kernel_s_ind,ksb,n,no_c);
    
    GminusI_ind = @(y) G_ind(y) - y;
    GTminusI_ind = @(y) GT_ind(y) - y;
    
    mid = sigma * GTminusI_ind(GminusI_ind(x1));
    mid1 = reshape(mid,[ksb n no_c]);
    
    p = x((abc-1)*ksb +1:abc*ksb,:,:);
    
    mid2 = (lambda*p) + adder;
    %mid2 = (1*p) + lambda*adder;
    mid_all = mid1+mid2;
    outp((abc-1)*ksb +1:abc*ksb,:,:) = mid_all;
end
outp = outp(:);
end


