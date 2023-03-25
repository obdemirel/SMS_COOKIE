function [oo] = atb_func(ddts,lls,zs,Ks,rho,Coil_sensitivites,ksb,slices)


[m,n,no_c] = size(ddts);

ET = @(x) encoder_t(x,m,n,no_c,slices,Coil_sensitivites,ksb);

zs = reshape(zs,m,n);
Ks = reshape(Ks,m,n);

for slis = 1:slices
ddts_all(:,:,:,slis) =  ddts(ksb*(slis-1) + 1:slis*ksb,:,:);
lambda_all(:,:,:,slis) =  lls(ksb*(slis-1) + 1:slis*ksb,:,:); 
z_all(:,:,:,slis) =  zs(ksb*(slis-1) + 1:slis*ksb,:,:);    
k_all(:,:,:,slis) =  Ks(ksb*(slis-1) + 1:slis*ksb,:,:); 
end


for slis = 1:slices
   insider_all(:,:,:,slis) = z_all(:,:,:,slis) - k_all(:,:,:,slis)/rho;
end

insiders = [];
for slis = 1:slices
   insiders = [insiders;insider_all(:,:,:,slis)];
end

insider = ET(insiders);

oo = [];
for slis = 1:slices
   oo = [oo;ddts_all(:,:,:,slis)+lambda_all(:,:,:,slis) + (rho/2)*insider(ksb*(slis-1) + 1:slis*ksb,:,:)];
end


oo = oo(:);

end

