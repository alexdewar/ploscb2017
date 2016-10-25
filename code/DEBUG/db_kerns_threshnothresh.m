clear

kname = 'vf_avkernels_r2';

load('vf_kernels_new.mat',kname)
kerns = eval(kname);

load('vf_kernels_nothresh.mat',kname)
kerns_nothresh = eval(kname);

clear(kname)

c = kstr_cents(kerns);
cth = kstr_cents(kerns_nothresh);

pcell = cell(1,2*numel(kerns));
for i = 1:numel(kerns)
    pcell{(i-1)*2+1} = [c(i,1), cth(i,1)];
    pcell{(i-1)*2+2} = [c(i,2), cth(i,2)];
end

% figure(1);clf
% plot(pcell{:})
% ksz = size(kerns(1).k);
% xlim([1 ksz(2)])
% ylim([1 ksz(1)])
% 
figure(3);clf
% subplot(1,2,1)
% showkernels(kerns(11),[],[],[-135 135],[-60 60])
% axis equal tight
% subplot(1,2,2)
showkernels(kerns_nothresh,[],[],[-135 135],[-60 60])
axis equal tight

% figure(1);clf
% for i = 1:numel(kerns)/2
%     subplot(2,7,i)
%     showkernels([kerns(i),kerns_nothresh(i)],[],.5)
% end

% figure(1);clf
% for i = 1:numel(kerns_nothresh)/2
%     subplot(2,7,i)
%     showkernels(kerns_nothresh(i),[],.5)
% end

% figure(2);clf
% for i = 1:numel(kerns_nothresh)
%     subplot(2,14,i)
%     showkernels(kerns(i),[],1)
% end