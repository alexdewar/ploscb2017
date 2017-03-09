clear
close all

% warning('off','images:imshow:magnificationMustBeFitForDockedFigure')

load('star001.mat');
% load('ofstad_etal_arena.mat')

tic
v = getviewfast_new(0,0,0,0,X,Y,Z,[],[]);
toc

% tic
% v2 = getView(0,0,0,0,X,Y,Z,false,false);
% toc

% v = true(10,100);
% x0 = 10;
% x1 = 50;
% y0 = 2;
% y1 = 8;
% 
% dx = x1-x0;
% dy = y1-y0;
% err = 0;
% if dx==0
%     v(y0:sign(dx):y1,x0) = false;
% else
%     derr = dy/dx;
%     y = y0;
%     for x = x0:sign(dx):x1
%         v(y,x) = false;
%         err = err+derr;
%         while err >= 0.5
%             v(y,x) = false;
%             y = y+sign(dx);
%             err = err-1;
%         end
%     end
% end

figure(2);clf
% subplot(2,1,1)
imshow(v)
% subplot(2,1,2)
% imshow(v2)
axis on