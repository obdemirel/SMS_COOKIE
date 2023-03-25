
clear all

    
load inital_ro_spirit.mat
dummy = new_kspace;
dummy(:,1:2:end,:,:) = 0;
dummy(:,77:100,:,:) = new_kspace(:,77:100,:,:);
new_kspace = dummy;

tic
kspace = new_kspace;
conc_acs_fregion = acs_f;
reg = 1e-2; %1e-2;

%%% KERNEL %%%%%%%%%
% XXX        XXX
% XXX    ->  X1X
% XXX        XXX
%%%%%%%%%%%%%%%%%%%%

kernel_row = 9;  %3
kernel_col = 9;  %3

kernel_dim = kernel_row*kernel_col;

MA = zeros((size(conc_acs_fregion,1)-(kernel_row-1))*(size(conc_acs_fregion,2)-(kernel_col-1)),kernel_dim*size(conc_acs_fregion,3));

%% MA matrix filling
for coil_selec = 1:size(conc_acs_fregion,3)
    selected_acs = conc_acs_fregion(:,:,coil_selec);
    row_count = 1;
        for col = 1:size(selected_acs,2)-(kernel_col-1)
            for row = 1:size(selected_acs,1)-(kernel_row-1)
                neighbors = selected_acs(row:row+(kernel_row-1),col:col+(kernel_col-1));
                neighbors = neighbors(:).';
                %neighbors_v = [neighbors(1:floor(kernel_dim/2)) 0 neighbors(ceil(kernel_dim/2)+1:end)];
                MA(row_count,(coil_selec-1)*(kernel_dim) +1:coil_selec*(kernel_dim)) = neighbors;
                row_count = row_count+1;
            end
        end
    disp(['MA ' num2str(coil_selec) ' coil is ready!'])
end


row_start = ceil(kernel_row/2);
row_end = size(conc_acs_fregion,1)-floor(kernel_row/2);
col_start = ceil(kernel_col/2);
col_end = size(conc_acs_fregion,2)-floor(kernel_col/2);


Mk = zeros(size(MA,1),size(conc_acs_fregion,3));
%% Mk vectors filling
for coil_selec = 1:size(conc_acs_fregion,3)
    selected_acs = conc_acs_fregion(row_start:row_end,col_start:col_end,coil_selec);
    Mk(:,coil_selec) = selected_acs(:);
end
disp('Mk is ready!')



lambda = norm(MA,'fro')/sqrt(size(MA,1))*reg;
ak = zeros(kernel_dim*size(conc_acs_fregion,3) -  1,size(conc_acs_fregion,3));

A = MA'*MA;
B = MA'*Mk;



for coil_selec = 1:size(conc_acs_fregion,3)
    tic
    
    A1 = A(1:ceil(kernel_dim/2)-1 + (coil_selec-1)*kernel_dim,:);
    A2 = A(ceil(kernel_dim/2)+1+ (coil_selec-1)*kernel_dim :end,:);
    newA = [A1;A2];
    A3 = newA(:,1:ceil(kernel_dim/2)-1 + (coil_selec-1)*kernel_dim);
    A4 = newA(:,ceil(kernel_dim/2)+1+ (coil_selec-1)*kernel_dim:end);
    newA = [A3 A4];
    
    B1 = B(1:ceil(kernel_dim/2)-1+ (coil_selec-1)*kernel_dim,:);
    B2 = B(ceil(kernel_dim/2)+1+(coil_selec-1)*kernel_dim:end,:);
    newB = [B1;B2];
    
    %ak(:,coil_selec) = newA\squeeze(newB(:,coil_selec));
    
    %ak(:,coil_selec) = inv(newA + eye(size(newA))*lambda)*squeeze(newB(:,coil_selec));
    ak(:,coil_selec) = (newA + eye(size(newA))*lambda)\squeeze(newB(:,coil_selec));
    
    disp(['Coil ' num2str(coil_selec) ' weights are ready!'])
    toc
end

new_ak = zeros(kernel_dim*size(conc_acs_fregion,3),size(conc_acs_fregion,3));


for coil_selec = 1:size(conc_acs_fregion,3)
    p1 = ak(1:ceil(kernel_dim/2) + (coil_selec-1)*kernel_dim -1,coil_selec);
    p2 = ak(ceil(kernel_dim/2) + (coil_selec-1)*kernel_dim:end,coil_selec);
    new_ak(:,coil_selec) = [p1;0;p2];
end

data_ak_set(:,:) = new_ak;
toc

display('Starting to reconstruction!')


for ss=1:size(new_kspace,4)
multi_slice_kspace = kspace(:,:,:,ss);
small_art_kspace1 = zeros(size(multi_slice_kspace,1)*slice_R,size(multi_slice_kspace,2),size(multi_slice_kspace,3),'single');
small_art_kspace1(1:slice_R:end,:,:) = multi_slice_kspace;
data_kspace(ss,:,:,:) = small_art_kspace1;
end

% data_ak(:,:) = ak;
kernel_s = kernel_row;
save('pre_data_golden','data_kspace','data_ak_set','kernel_s','-v7.3')
