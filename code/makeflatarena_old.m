
ht = 1; %128e-3;
d = 123e-3;
wd_bar = ht/3;
nvert = 4;

dmatsz = 1000;

% wd = pi*d;
wd = 1;

% vertical bars
Xvert = repmat([0 0 1 1]*wd_bar,nvert,1);
Yvert = repmat([0 1 1 0]*ht,nvert,1);
Xvert = Xvert+repmat(2*wd_bar*(0:nvert-1)',1,4);

% horizontal bar
Xhorz = [0 1 1 0]*wd/3;
Yhorz = [1 1 2 2]*wd_bar;

% 45 deg bars
k = ht/wd;
chg = wd_bar/(2*sqrt(2));
% xchg = wd_bar/sqrt(k^2+1);
% ychg = k*xchg;
ychg = wd_bar*sin(pi/4);
ystart = 0.1; % ht-ychg;
xstart = 0;

dmat = zeros(dmatsz);
dbarwd = dmatsz*wd_bar/2;
[xi,yi] = meshgrid(0:dmatsz-1,0:dmatsz-1);
% figure(1)
c = 3*sqrt(2)*dbarwd*[-1 0 1];
for i = 1:length(c)
    dmat(abs(xi-yi+c(i))<=dbarwd) = 1;
%     imshow(dmat)
end
figure(1);clf
[X,Y] = bw2polygon(dmat);
alfill(X,Y,'k');


% x45 = xstart+[0 -xchg/2 -xchg wd-xchg wd-xchg/2 wd];
% y45 = ystart+[0 ychg/2 ychg ht+ychg ht-ychg/2 ht];
% x45 = xstart+[0 -chg/2 -chg wd wd+chg/2 wd+chg];
% y45 = ystart+[0 chg/2 chg ht+2*chg ht+1.5*chg ht+chg];
% x45 = xstart+[0 -chg 0 wd-wd_bar/sin(pi/4) wd wd+chg];
% y45 = ystart+[0 chg 2*chg ht ht+2*ychg ht+ychg];
% x45 = xstart+[0 0 wd-wd_bar/sin(pi/4) wd];
% y45 = ystart+[0 chg ht ht];

% x45 = [0 -chg wd wd+chg]';
% y45 = [0 chg ht+2*chg ht+chg]';
% hts = 0:2*ychg:ht;
% hts = [hts -hts(2:end)];
% y45 = bsxfun(@plus,y45,hts);
% % y45 = [bars,y45,-bars];
% x45 = repmat(x45,1,size(y45,2));
% [x45m,y45m] = boxshape(x45,y45);

% x45(1) = 0;
% y45(2) = ystart;
% yout = wd+ystart;
% if yout > ht
%     y45(2) = ht;
%     x45(2) = ht-ystart;
% else
%     y45(2) = yout;
%     x45(2) = wd;
% end
% y45(3) = yout + 2*chg;
% if y45(3) > ht
%     y45(3) = ht;
%     x45(3) = ht-ystart-2*chg;
% else
%     x45(3) = wd;
% end
% y45(4) = ystart+2*chg;
% if y45(4) > ht
%     y45(4) = ht;
%     x45(4) = ht-ystart-2*chg;
% else
%     x45(4) = 0;
% end


% x45m = max(0,min(wd,x45));
% y45m = max(0,min(ht,y45));

% figure(1);clf
% fill(x45,y45,'r');
% axis equal
% 
% figure(2);clf
% fill(x45m,y45m,'k')
% axis equal