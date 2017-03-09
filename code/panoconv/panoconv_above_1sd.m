function [acts,angs]=panoconv_above_1sd(im,rkerns,fov,nangs)
if nargin < 3;
    fov = 360;
end
if nargin < 4
    nangs = size(im,2);
end

% how much to trim from image
xoff = round(0.5*size(im,2)*(1-fov/360));
xrng = size(im,2);
xvals = round(linspace(xrng/2,-xrng/2-1,nangs));

% if nangs==size(im,2)
%     
% else
    acts = NaN(nangs,1);
    for i = 1:nangs
        rim = circshift(im,[0 xvals(i)]);
        cim = rim(:,1+xoff:end-xoff);

        cacts = getacts(cim,rkerns);
        acts(i) = mean(cacts > std(cacts));
    end
% end

if nargout==2
    angs = linspace(-180,180,nangs);
end