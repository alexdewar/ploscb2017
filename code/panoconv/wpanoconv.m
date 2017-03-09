function [acts,angs,act0]=wpanoconv(im,rkerns,weights,fov,nangs)
if nargin < 4;
    fov = 360;
end
if nargin < 5
    nangs = size(im,2);
end

weights = weights(:);

% how much to trim from image
xoff = round(0.5*size(im,2)*(1-fov/360));
xrng = size(im,2);
xvals = round(linspace(xrng/2,-xrng/2-1,nangs));

angs = linspace(-180,180,nangs+1);
angs = angs(1:end-1);
sweights = sum(weights);
acts = NaN(nangs,1);
% [~,i90] = min(abs(angs+90));
for i = 1:nangs
    rim = circshift(im,[0 xvals(i)]);
    cim = rim(:,1+xoff:end-xoff);

    cacts = weights.*(1+getacts(cim,rkerns))/2;
    acts(i) = sum(cacts)/sweights;
    if xvals(i)==0
        act0 = cacts;
        
%         figure(10);clf
%         subplot(3,1,1)
%         imshow(im)
%         subplot(3,1,2)
%         imshow(cim)
%         subplot(3,1,3)
%         imshow(im90)
%         keyboard
%     elseif i==i90
%         act90 = cacts;
        
%         im90 = cim;
    end
end