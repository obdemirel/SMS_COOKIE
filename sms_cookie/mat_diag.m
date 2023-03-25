function m = mat_diag(x,m,n,no_c,data_ak_ind,kernel_s_ind,ksb,slices,sigmas)


x = reshape(x,m,n,no_c);

for abc = 1:slices
    ak_ind = squeeze(data_ak_ind(abc,:,:));
    sigma = sigmas(abc);
    
    x1 = x((abc-1)*ksb +1:abc*ksb,:,:);
    x1 = x1(:);
    
    G_ind = @(x1) conv_op_ind(ak_ind,x1,kernel_s_ind,kernel_s_ind,ksb,n,no_c);
    GT_ind = @(x1) conv_op_ind_t(ak_ind,x1,kernel_s_ind,kernel_s_ind,ksb,n,no_c);
    
    GminusI_ind = @(x1) G_ind(x1) - x1;
    GTminusI_ind = @(x1) GT_ind(x1) - x1;
    
    mid = sigma * GTminusI_ind(GminusI_ind(x1));
    k((abc-1)*ksb +1:abc*ksb,:,:) = reshape(mid,ksb,n,no_c);
    
end
m = k(:);
end

