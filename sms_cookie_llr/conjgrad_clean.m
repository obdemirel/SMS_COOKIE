function [x] = conjgrad_clean(niter,ATA,b, x)
r = b - ATA(x);
p = r;
rsold = r' * r;
% figure,
for i = 1:length(b)
    Ap = ATA(p);
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;
    if sqrt(rsnew) < 1e-25 
        disp(['Ended because error too small!'])    
        break;
    end
    if i>niter
        disp(['Ended because iteration ends!'])    
        break;
    end
    p = r + (rsnew / rsold) * p;
%     error(i) = norm(ATA(x)-b);
    rsold = rsnew;
    
    
    %disp(['iteration: ' num2str(i) ' and error: ' num2str(error(i))])
    %%%%%%%%%%%PLOT%%%%%%%%%%%%%%%%%%%%%%%
%     i
%      rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);
%      asd = reshape(x,[480 152 30]);
%      asd = flipud(rssq(ifft2(asd(161:320,:,:))));
%      asd = asd./max(asd(:));
%      subplot(1,3,1), imshow(asd,[0 1]), drawnow()
%      subplot(1,3,2), plot((error)), drawnow()
%      subplot(1,3,3), plot(log(error)), drawnow()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
end