function irf = makeIRF(pixperdeg, fov, ...
                       area_ex,  ecc_ex, cent_ex, th_ex, ...
                       propr_in, ecc_in, cent_in, th_in, ...
                       propr_0,  ecc_0,  cent_0,  th_0)

if isempty(pixperdeg)
    pixperdeg = 1;
end

if isempty(fov)
    fov = [-60 60; -135 135];
end

fprintf(['area_ex %f\necc_ex %f\ncent_ex [%f %f]\nth_ex %f\n' ...
         'propr_in %f\necc_in %f\ncent_in [%f %f]\nth_in %f\n' ...
         'propr_0 %f\necc_0 %f\ncent_0 [%f %f]\nth_0 %f\n\n'], ...
          area_ex,  ecc_ex, cent_ex(1), cent_ex(2), th_ex, ...
          propr_in, ecc_in, cent_in(1), cent_in(2), th_in, ...
          propr_0,  ecc_0,  cent_0(1), cent_0(2),  th_0);

% % excitatory
% area_ex = 1000; % deg^2
% ecc_ex = 0.5;
% cent_ex(1) = 0; % deg
% cent_ex(2) = 20; % deg
% th_ex = pi/2; % rad
% 
% % inhibitory
% propr_in = 2;
% ecc_in = 0.5;
% cent_in(1) = 0;
% cent_in(2) = 0;
% th_in = pi/2;
% 
% % null
% propr_0 = 1.6;
% ecc_0 = 0.5;
% cent_0(1) = 0;
% cent_0(2) = 0;
% th_0 = 0;

im_size = round(pixperdeg*range(fov,2)');

% get image xs & ys (deg)
im = zeros(im_size);
% [imx,imy] = ndgrid((1:im_size(2)) + cent_ex(2), ...
%                    (1:im_size(1)) + cent_ex(1));
[imx,imy] = meshgrid(linspace(fov(2,1),fov(2,2),im_size(2))-cent_ex(1), ...
                     linspace(fov(1,2),fov(1,1),im_size(1)));
% [imx,imy] = rotatexy(imx,imy,-th_ex);
% imx = imx + cent_ex(2);
% imy = imy + cent_ex(1);

% inhibitory region
% sni = sin(th_in);
% csi = cos(th_in);
% ki = sqrt(1-ecc_in.^2);
% im( (imx*csi - imy*sni - cent_in(2)).^2 * ki + ...
%     (imx*sni + imy*csi - cent_in(1)).^2 / ki ...
%                             <= propr_in.^2 * area_ex/pi ) = -1;
% 
% % null region
% sn0 = sin(th_0);
% cs0 = cos(th_0);
% k0 = sqrt(1-ecc_0.^2);
% im( (imx*cs0 - imy*sn0 - cent_0(2)).^2 * k0 + ...
%     (imx*sn0 + imy*cs0 - cent_0(1)).^2 / k0 ...
%                             <= propr_0.^2 * area_ex/pi ) = 0;

% excitatory region
ke = sqrt(1-ecc_ex.^2);
pos = (imx.^2 * ke + imy.^2 / ke) <= (area_ex/pi);
im(pos) = 1/sum(pos(:));

% normalise inhibitory vals
neg = im < 0;
im(neg) = -1/sum(neg(:));

foff = [fov(2,1),fov(1,1)];
irf = struct('k',im,'cent',pixperdeg*(cent_ex-foff), ...
             'ex', struct('area',pixperdeg.^2*area_ex,'ecc',ecc_ex,'th',th_ex), ...
             'in', struct('propr',propr_in,'ecc',ecc_in,'cent',pixperdeg*(cent_in-foff),'th',th_in), ...
             'null', struct('propr',propr_0,'ecc',ecc_0,'cent',pixperdeg*(cent_0-foff),'th',th_0));
        
% % show
% figure(1);clf
% showkernel(irf)
% set(gca,'YDir','normal');