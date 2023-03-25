function y = a_operator_t(x,m, n, no_c,ksb)

x = reshape(x,m,n,no_c);


% for co =1:no_c
% im_all(:,:,co) = ifft2(x(:,:,co));
% end

im_all = ifft2(x);


% for co =1:no_c
%     for slice_number = 1:3
%         y((slice_number-1)*ksb +1:slice_number*ksb,:,:) = fft2(im_all((slice_number-1)*ksb +1:slice_number*ksb,:,:));
%     end
% end

% for co =1:no_c
    for slice_number = 1:3
        y((slice_number-1)*ksb +1:slice_number*ksb,:,:) = fft2(im_all((slice_number-1)*ksb +1:slice_number*ksb,:,:));
    end
% end


y = y(:);

end