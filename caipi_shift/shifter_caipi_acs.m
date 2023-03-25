function [acs] = shifter_caipi_acs(acs,fov_shifts,fov_amount);

for ii=1:size(acs,4)
[m,n,no_c] = size(squeeze(acs(:,:,:,ii)));
phases = exp(sqrt(-1)*(fov_shifts(ii))*2*pi/fov_amount*(0:n-1));
phas = repmat(phases,[m 1 no_c]);
acs(:,:,:,ii) = squeeze(acs(:,:,:,ii)) .* phas;
end

end

