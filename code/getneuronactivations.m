function [act,tkerns]=getneuronactivations(ims,tkerns)
if nargin < 2
    tkerns = 'both';
end
    
if ischar(tkerns)
    if strcmpi(tkerns,'both')
        load('vf_tkernels.mat','tkerns_r2','tkerns_r4');
        tkerns = [tkerns_r2,tkerns_r4];
    else
        vname = ['tkerns_' tkerns];
        load('vf_tkernels.mat',vname);
        tkerns = eval(vname);
    end
end

vsz = [size(ims,1),size(ims,2)];

act = NaN(numel(tkerns),1,size(ims,3));
for i = 1:numel(tkerns)
    if ~all(vsz./tkerns(i).im_size == 1)
        tkerns(i) = resizekernel(tkerns(i),vsz,0.25);
    end
    [skern,yind,xind] = shiftkern(tkerns(i),vsz,tkerns(i).oldcent);
    act(i,1,:) = shiftdim(sum(sum(bsxfun(@times,skern,ims(yind,xind,:)))));
    
%     figure(1);clf
%     showkernel(skern)
%     drawnow
%     pause(0.5);
end