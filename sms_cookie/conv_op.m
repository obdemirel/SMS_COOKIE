function [new_kspace] = conv_op(ak,mb_kspace,kernel_r,kernel_c,m,n,coil_n,acq_p)
    mb_kspace = reshape(mb_kspace,m,n,coil_n);
    new_kspace = zeros(size(mb_kspace),'single');
    for coil_sel = 1:size(mb_kspace,3)
        %tic
        selected_coil_ak = ak(:,coil_sel);
        selected_coil_ker = reshape(selected_coil_ak,kernel_r,kernel_c,size(mb_kspace,3));
        convelved_space = zeros(size(mb_kspace,1),size(mb_kspace,2));
        for coil_conv = 1:size(mb_kspace,3)
            selected_coil = ((mb_kspace(:,:,coil_conv)));
            selected_ker =  selected_coil_ker(:,:,coil_conv);
            kernel_2D = selected_ker;
            c_space = filter2(kernel_2D,selected_coil);
            %c_space = conv2(selected_coil,kernel_2D,'same');
            convelved_space = convelved_space + c_space;
        end
        new_kspace(:,:,coil_sel) = (convelved_space);
        %display(['Recon Coil=' num2str(coil_sel) ' is ready'])
        %toc
    end
    %new_kspace(acq_p) = mb_kspace(acq_p);
    new_kspace = new_kspace(:);
end

