%tic
clear all
%close all

for klm=1:5
load inital
load inital_ind


kspace = new_kspace;
conc_acs_fregion = ACS1_final+ACS2_final+ACS3_final;
reg = 1e-2;

%%% KERNEL %%%%%%%%%
% XXX        XXX
% XXX    ->  X1X
% XXX        XXX
%%%%%%%%%%%%%%%%%%%%

kernel_row = 7;  %3
kernel_col = 7;  %3

kernel_dim = kernel_row*kernel_col;

MA = zeros((size(conc_acs_fregion,1)-(kernel_row-1))*(size(conc_acs_fregion,2)-(kernel_col-1)),kernel_dim*size(kspace,3));

%% MA matrix filling
for coil_selec = 1:size(conc_acs_fregion,3)
    selected_acs = conc_acs_fregion(:,:,coil_selec);
    row_count = 1;
        for col = 1:size(selected_acs,2)-(kernel_col-1)
            for row = 1:size(selected_acs,1)-(kernel_row-1)
                neighbors = selected_acs(row:row+(kernel_row-1),col:col+(kernel_col-1));
                neighbors = neighbors(:).';
                MA(row_count,(coil_selec-1)*(kernel_dim) +1:coil_selec*(kernel_dim)) = neighbors;
                row_count = row_count+1;
            end
        end
    disp(['MA ' num2str(coil_selec) ' coil is ready!'])
end

if(klm==1)
    conc_acs_fregion = ACS1_final;
end
if(klm==2)
    conc_acs_fregion = ACS2_final;
end
if(klm==3)
    conc_acs_fregion = ACS3_final;
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


reg = norm(MA,'fro')/sqrt(size(MA,1))*reg;
ak = zeros(kernel_dim*size(kspace,3),size(kspace,3));
 
for coil_selec = 1:size(conc_acs_fregion,3)
    tic
	MA(:,floor(kernel_dim/2)+((coil_selec-1)*kernel_dim)) = 0;
    %ak(:,coil_selec) = pinv(MA'*MA)*MA'*Mk(:,coil_selec);
    A = MA'*MA + reg*eye(size(MA,2));
    b = MA'*Mk(:,coil_selec);
    ak(:,coil_selec) = A\b;
    disp(['Coil ' num2str(coil_selec) ' weights are ready!'])
    toc
end

clear MA
clear Mk

data_ak_set(klm,:,:) = ak;
end

display('Starting to reconstruction!')


for ss=1:40
multi_slice_kspace = kspace(:,:,:,ss);
small_art_kspace1 = zeros(size(multi_slice_kspace,1)*slice_R,size(multi_slice_kspace,2),size(multi_slice_kspace,3));
small_art_kspace1(1:slice_R:end,:,:) = multi_slice_kspace;
data_kspace(ss,:,:,:) = small_art_kspace1;
end

% data_ak(:,:) = ak;
kernel_s = kernel_row;
save('pre_data','data_kspace','data_ak_set','kernel_s','-v7.3')


% save('kspace_data','xk_prev')
% save('orig_kspace','small_art_kspace1');
% save('spirit_kernels','ak')
% save('acs_kspace','conc_acs_fregion')
% save('kernel_s','kernel_row')