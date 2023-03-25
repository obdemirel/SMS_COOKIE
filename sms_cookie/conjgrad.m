function [x,error] = conjgrad(niter,ATA,b, x,gui_on)
r = b - ATA(x);
p = r;
rsold = r' * r;

for i = 1:length(b)
    if(gui_on==1)
    progressbar(i/(niter+1))
    end
    Ap = ATA(p);
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;
    if sqrt(rsnew) < 1e-15
        disp(['Ended because error too small!'])
        break;
    end
    if i>niter
        disp(['Ended because iteration ends!'])
	%error(i) = 1;
	%x
        break;
    end
    p = r + (rsnew / rsold) * p;
    %%% remove the comment to calculate the residual error
    %%% not it increases the reconstruction time
    error(i) = 1;%norm(b-ATA(x))/norm(b);
    rsold = rsnew;
    
%     asd = reshape(x,[480 144 30]);
%     rssq = @(x) squeeze(sum(abs(x).^2,3)).^(1/2);
%     figure, imshow(rssq(ifft2c(asd(161:2*160,:,:))),[])
%     figure, imshow(log((asd(161:2*160,:,10))),[])
end
if(gui_on==1)
progressbar(1)
end
end