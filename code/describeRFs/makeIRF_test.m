% function irf = makeIRF_test(im_size, ...
%                        area_ex,  ecc_ex, cent_ex, th_ex, ...
%                        propr_in, ecc_in, cent_in, th_in, ...
%                        propr_0,  ecc_0,  cent_0,  th_0)
% 
% fprintf(['area_ex %f\necc_ex %f\ncent_ex [%f %f]\nth_ex %f\n' ...
%          'propr_in %f\necc_in %f\ncent_in [%f %f]\nth_in %f\n' ...
%          'propr_0 %f\necc_0 %f\ncent_0 [%f %f]\nth_0 %f\n\n'], ...
%           area_ex,  ecc_ex, cent_ex(1), cent_ex(2), th_ex, ...
%           propr_in, ecc_in, cent_in(1), cent_in(2), th_in, ...
%           propr_0,  ecc_0,  cent_0(1), cent_0(2),  th_0);

im_size = [120 270];

% excitatory
area_ex = 250; % deg^2
ecc_ex = 0.5;
cent_ex(1) = 60; % deg
cent_ex(2) = 60; % deg
th_ex = pi/2; % rad

% inhibitory
propr_in = 2;
ecc_in = 0.5;
cent_in(1) = 0;
cent_in(2) = 0;
th_in = pi/2;

% null
propr_0 = 1;
ecc_0 = 0.5;
cent_0(1) = 0;
cent_0(2) = 0;
th_0 = 0;

% get image xs & ys (deg)
im = zeros(im_size);
[imx,imy] = meshgrid((1:im_size(2))-cent_ex(1), (1:im_size(1)) - cent_ex(2));
[imx,imy] = rotatexy(imx,imy,th_ex);

% inhibitory region
sni = sin(th_in);
csi = cos(th_in);
ki = sqrt(1-ecc_in.^2);
sn0 = sin(th_0);
cs0 = cos(th_0);
k0 = sqrt(1-ecc_0.^2);
inh = ((1/(propr_in.^2 * area_ex))*((imx*csi - imy*sni - cent_in(2)).^2 * ki + ...
      (imx*sni + imy*csi - cent_in(1)).^2 / ki)) - ...
      ((1/(propr_0.^2 * area_ex))*((imx*cs0 - imy*sn0 - cent_0(2)).^2 * k0 + ...
      (imx*sn0 + imy*cs0 - cent_0(1)).^2 / k0));
      
im(inh==0) = -1;

figure(10);clf
surf(inh)

% excitatory region
ke = sqrt(1-ecc_ex.^2);
pos = (imx.^2 * ke + imy.^2 / ke) <= (area_ex/pi);
im(pos) = 1/sum(pos(:));

% normalise inhibitory vals
neg = im < 0;
im(neg) = -1/sum(neg(:));

irf = struct('k',im,'cent',cent_ex, ...
             'ex', struct('area',area_ex,'ecc',ecc_ex,'th',th_ex), ...
             'in', struct('propr',propr_in,'ecc',ecc_in,'cent',cent_in,'th',th_in), ...
             'null', struct('propr',propr_0,'ecc',ecc_0,'cent',cent_0,'th',th_0));
        
% show
figure(1);clf
showkernel(irf)