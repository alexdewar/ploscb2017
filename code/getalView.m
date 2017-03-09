function v=getalView(x,y,z,th,X,Y,Z,av,dohisteq,doclf)
%function v=getView(x,y,z,th,X,Y,Z,av,dohisteq)

if nargin < 10
    doclf = false;
end
if nargin < 9
    dohisteq = false;
end
if nargin < 8
    av = true;
end

% get the ground height
z0=0;%getHeight(x,y,X,Y,Z);

ff=figure(999);
if doclf
    clf
end

set(ff,'Position',[3,524,1044,304])
% If you need to resize the figure that pops up reset this line
% by first resizing the image and then typing get(gcf,'Position') to return
% the current position. Put those values into the line above.
[TH,PHI,R]=cart2sph(X-x,Y-y,Z-z-z0);

% remove small far objects for speed up
% [t,r]=cart2pol(mean(R,2),mean(Z,2));
% TH=TH(t/2/pi*360>1,:);
% PHI=PHI(t/2/pi*360>1,:);

TH=pi2pi(TH-th);

ind=(max(TH')-min(TH')<pi);
alfill(-TH(ind,:)',PHI(ind,:)','k')
hold on

negind=(max(TH')-min(TH')>pi);
T2=TH(negind,:);
P2=PHI(negind,:);

TP=T2;
TP(T2<0)=TP(T2<0)+2*pi;

alfill(-TP',P2','k')

TP=T2;
TP(T2>0)=TP(T2>0)-2*pi;

alfill(-TP',P2','k')

set(gca,'cameraViewAngle',1)

axis equal
axis([-pi pi 0 1.2])
axis off
hold off
drawnow

f=getframe;
% this function needs the image processing toolbox
v=imresize(f.cdata,900/size(f.cdata,2));
v=double(v(:,:,1)>0);
%v=flipud(v);

if av
    % this function needs more recent versions of MATLAB
    v=blockproc(v,[10 10],@(x)mean2(x.data));
    v=v(1:end-1,:);
else
    v=v(1:end-4,:);
end

if dohisteq
    v = histeq(v);
end

% function m=mean1(X)
% m=mean(X(:));

function x=pi2pi(x)
x=mod(x,2*pi);
x=x-(x>pi)*2*pi;