function [x_coils_alt] = LOST4(x_coils,threshold_factor,indices)
addpath('Utilities')

%% GET THE SIZE AND UNDERSAMPLING INFO
no_c = 1;

%% LEARN SIMILARITY CLUSTERS
Nb = 4; Ngroup = 16; Nsearch = 8; Ndepth = 1;
fprintf('\n'); fprintf('Running LOST stage 2'); fprintf('\n');
match_threshold = 0.05; % this can be changed between 0-1.
% for ind1 = 1:no_c
%     tic
%     fprintf('Finding similarity clusters in coil %d', ind1);
%     fprintf('... This may take a while...'); fprintf('\n');
%     mm = x_coils(:,:,:,ind1);
%     im_train = x_coils(:,:,:,ind1)/max(abs(mm(:))); im_train = abs(im_train);
%     GroupIndices_cbc{ind1} = FindSimilarBlocksInC3D(single(im_train), single(match_threshold), single(Nb), single(Nsearch), single(1), single(Ngroup), single(Ndepth));
%     toc
% end

GroupIndices_cbc = indices;

beta_kaiser = 2; w_kaiser2d = kaiser(Nb, beta_kaiser) * kaiser(Nb,beta_kaiser)';

for ind2 = 1:no_c
    x =  x_coils(:,:,:,ind2);
    URange = max(abs(x(:)));
    GroupIndices =  GroupIndices_cbc{ind2};
    for ind1 = 1:15
        if mod(ind1, 2) == 0,
            optstruct.scheme = single(0); optstruct.thresh = single(threshold_factor * URange);
        else
            optstruct.scheme = single(2); optstruct.thresh = single(threshold_factor * URange)^2;
        end
        C_Output = BlockThresholdingInC3D(single(x), optstruct, GroupIndices, single(Nb), single(w_kaiser2d), single(1));
        x = C_Output;
    end
    x_coils_alt(:,:,:,ind2) = x;    %This is data-consistent
end

end

