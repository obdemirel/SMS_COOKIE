%tic
clear all
%close all


load('C:\Users\obura\Desktop\dMRI\spsg\initial\inital_ind.mat')



for abc = 1:4
    if abc ==1
        conc_acs_fregion = ACS1_final;
    elseif abc == 2
        conc_acs_fregion = ACS2_final;
    elseif abc == 3
        conc_acs_fregion = ACS3_final;
        elseif abc == 4
        conc_acs_fregion = ACS4_final;
    end
    
    
    reg = 1e-2;
    
    %%% KERNEL %%%%%%%%%
    % XXX        XXX
    % XXX    ->  X1X
    % XXX        XXX
    %%%%%%%%%%%%%%%%%%%%
    
    kernel_row = 5;  %3
    kernel_col = 5;  %3
    
    kernel_dim = kernel_row*kernel_col;
    
    MA = zeros((size(conc_acs_fregion,1)-(kernel_row-1))*(size(conc_acs_fregion,2)-(kernel_col-1)),(kernel_dim-1)*size(conc_acs_fregion,3));
    
    %% MA matrix filling
    for coil_selec = 1:size(conc_acs_fregion,3)
        selected_acs = conc_acs_fregion(:,:,coil_selec);
        row_count = 1;
        for col = 1:size(selected_acs,2)-(kernel_col-1)
            for row = 1:size(selected_acs,1)-(kernel_row-1)
                neighbors = selected_acs(row:row+(kernel_row-1),col:col+(kernel_col-1));
                neighbors = neighbors(:).';
                neighbors_v = [neighbors(1:floor(kernel_dim/2)) neighbors(ceil(kernel_dim/2)+1:end)];
                MA(row_count,(coil_selec-1)*(kernel_dim-1) +1:coil_selec*(kernel_dim-1)) = neighbors_v;
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
    
    
    
    reg = norm(MA,'fro')/sqrt(size(MA,1))*reg;
    ak = zeros((kernel_dim-1)*size(conc_acs_fregion,3),size(conc_acs_fregion,3),'single');
    
    
    for coil_selec = 1:size(conc_acs_fregion,3)
        tic
        %     MA(:,floor(kernel_dim/2)+((coil_selec-1)*kernel_dim)) = 0;
        %ak(:,coil_selec) = pinv(MA'*MA)*MA'*Mk(:,coil_selec);
        A = MA'*MA + reg*eye(size(MA,2));
        b = MA'*Mk(:,coil_selec);
        ak(:,coil_selec) = A\b;
        disp(['Coil ' num2str(coil_selec) ' weights are ready!'])
        toc
    end
    
    new_ak = zeros((kernel_dim)*size(conc_acs_fregion,3),size(conc_acs_fregion,3),'single');

    cons1 = kernel_dim-1;
    cons2 = kernel_dim;
    for coil_selec = 1:size(conc_acs_fregion,3)
        dumm1 = ak((coil_selec-1)*cons1 +1:coil_selec*cons1,:);
        dumm2 = [dumm1(1:cons1/2,:);zeros(1,size(conc_acs_fregion,3));dumm1(cons1/2 + 1:end,:)];
        new_ak((coil_selec-1)*cons2 +1:coil_selec*cons2,:) = dumm2;
    end
    
    clear MA
    clear Mk
    
    display('Starting to reconstruction!')
    
    ak = new_ak;
    
    data_ak_ind(abc,:,:) = ak;
    kernel_s_ind = kernel_row;
    
end


save('pre_data_ind_fixed2','data_ak_ind','kernel_s_ind')

% save('kspace_data','xk_prev')
% save('orig_kspace','small_art_kspace1');
% save('spirit_kernels','ak')
% save('acs_kspace','conc_acs_fregion')
% save('kernel_s','kernel_row')