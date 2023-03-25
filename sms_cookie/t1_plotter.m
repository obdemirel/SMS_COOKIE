
clear all
load data_new_t1
%%% for spirit phase sense

figure,
suptitle('SPIRIT Only')

t1Map_spirit_phase_sens(:,:,2) = circshift(t1Map_spirit_phase_sens(:,:,2),-50,2);
t1Map_spirit_phase_sens(:,:,3) = circshift(t1Map_spirit_phase_sens(:,:,3),-100,2);

for i=1:3

subplot(1,3,i), imagesc((t1Map_spirit_phase_sens(40:110,25:95,i)), [0 3000]), colormap(ckmapLLCine), colorbar
%subplot(1,3,i), imagesc((t1Map_spirit_phase_sens(:,:,i)), [0 3000]), colormap(ckmapLLCine), colorbar
axis square;
if(i==1)
    ylabel('Phase Sense')
end
end



load data_new_t1_sg
figure,
suptitle('SENSE-GRAPPA')

t1Map_grappa_phase_sens(:,:,2) = circshift(t1Map_grappa_phase_sens(:,:,2),-50,2);
t1Map_grappa_phase_sens(:,:,3) = circshift(t1Map_grappa_phase_sens(:,:,3),-100,2);

for i=1:3

subplot(1,3,i), imagesc((t1Map_grappa_phase_sens(40:110,25:95,i)), [0 3000]), colormap(ckmapLLCine), colorbar
%subplot(1,3,i), imagesc((t1Map_grappa_phase_sens(:,:,i)), [0 3000]), colormap(ckmapLLCine), colorbar
axis square;
if(i==1)
    ylabel('Phase Sense')
end
end



load reg_t1
figure,
suptitle('SENSE-GRAPPA')

t1Map_spirit_phase_sens(:,:,2) = circshift(t1Map_spirit_phase_sens(:,:,2),-50,2);
t1Map_spirit_phase_sens(:,:,3) = circshift(t1Map_spirit_phase_sens(:,:,3),-100,2);

for i=1:3

subplot(1,3,i), imagesc((t1Map_spirit_phase_sens(40:110,25:95,i)), [0 3000]), colormap(ckmapLLCine), colorbar
%subplot(1,3,i), imagesc((t1Map_grappa_phase_sens(:,:,i)), [0 3000]), colormap(ckmapLLCine), colorbar
axis square;
if(i==1)
    ylabel('Phase Sense')
end
end