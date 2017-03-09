function vf_panoconv(im,r2_or_r4,glomnum_for_averagekernel,use_left_version_of_kernel,save_figures)
    if nargin < 5
        save_figures = true;
    end
    if nargin < 4
        use_left_version_of_kernel = true;
    end
    if nargin < 3
        glomnum_for_averagekernel = 5;
    end
    if nargin < 2
        r2_or_r4 = 'r4';
    end
    if nargin < 1
        im = 'dartboard.mat';
    end
    

    if ischar(im)
        im = im2double(loadim((im)));
    end

%     r2_or_r4 = 'r2';
%     glomnum_for_averagekernel = 1;
%     use_left_version_of_kernel = true;
%     image_filename = 'dartboard.mat';
%     save_figures = true;

    thresh = 0.25;
    krng = [-60 60; -135 135]';
    irng = [-52 52; -180 180]';

    load('vf_kernels.mat');
%     kernels = eval(['vf_kernels_' r2_or_r4]);
    avkernels = eval(['vf_avkernels_' r2_or_r4]);
    whk = celleq(glomnum_for_averagekernel,avkernels.glomnum) & ...
          celleq(use_left_version_of_kernel,avkernels.isleft);
    kernel = avkernels(whk);
    % kernels to use to calculate where centre should be for average kernel
%     compkernels = kernels(celleq(glomnum_for_averagekernel,kernels.glomnum));
    
    %% resize kernel
%     im = im2double(rgb2gray(imread(image_filename)));
%     load(image_filename);
%     im = im2double(dartboardim);
    kern = (sign(kernel.k)+1)/2; % set kernel to be between 0 and 1

    % calculate new size
    ifov = range(irng);
    kfov = [ifov(1),range(krng(:,2))]; % squish kernel vertically
    origksize = size(kern);
    szfac = size(kern).*(kfov.*size(im))./(origksize.*ifov);
    kern = imresize(kern,szfac); % resize kernel
    
    kern = (kern*2)-1; % set kernel to be between -1 and 1
    kern(abs(kern)<thresh) = 0; % remove sub-threshold values
    kern(kern>0) = 1/sum(sum(kern>0)); % normalise kernel values
    kern(kern<0) = -1/sum(sum(kern<0));

    % figure(1);clf
    % imagesc(sign(kernel));
    % colormap(neuroncolormap);

    midi = floor(size(kern,2)/2);
    oldimwd = size(im,2);
    vals = NaN(oldimwd,1);
    im2 = [im(:,end-midi+1:end) im im(:,1:size(kern,2)-midi)];

    % figure(1);clf
    % imshow(im)
    % return
%     avkcent = scalecent(kernel.cent,origksize,size(kern));
%     kcents = cell2mat(cellfun(@transpose,{compkernels.cent},'uniformoutput',false))';
%     kcents = scalecent(kcents,repmat(origksize,size(kcents,1)), ...
%                               repmat(size(kern),size(kcents,1)));
    
    ths = linspace(irng(1,2),irng(2,2),length(vals));
    ymax = 0;
    ymin = inf;
%     startprogbar(20,numel(compkernels)*length(vals),'Doing convolutions...');
%     for i = 1:numel(compkernels)
%         coff = round(kcents(i,:)-avkcent);
%         if coff(2) < 0
%             ckern = [kern(1-coff(2):end,:);zeros(-coff(2),size(kern,2))];
%         else
%             ckern = [zeros(coff(2),size(kern,2));kern(1:end-coff(2),:)];
%         end
%         if coff(1) < 0
%             ckern = [ckern(:,1-coff(1):end) zeros(size(ckern,1),-coff(1))];
%         else
%             ckern = [zeros(size(ckern,1),coff(1)) ckern(:,1:end-coff(1))];
%         end

%         figure(1);clf
%         subplot(1,3,1);
%         showkernel(imresize(compkernels(i).k,szfac),kcents(i,:));
%         subplot(1,3,2);
%         showkernel(kern,avkcent);
%         subplot(1,3,3);
%         showkernel(ckern,kcents(i,:));
%         keyboard

        for j = 1:length(vals)
            vals(j) = sum(sum(kern.*im2(:,j-1+(1:size(kern,2)))));
            
%             progbar; %(sprintf('Using kernel %d/%d... ',i,numel(compkernels)));
        end
        
        maxval = max(vals);
        if maxval > ymax
            ymax = maxval;
        end
        minval = min(vals);
        if minval < ymin
            ymin = minval;
        end
        
        figure(1);clf
        plot(ths,vals);
        axis tight

%         ksz = size(compkernels(i).k)-1;
%         prop = (kcents(i,[2 1])-1)./ksz;
%         ccoord = prop.*ifov+irng(1,:);
%         title(sprintf('Centre at (%.2fdeg, %.2fdeg); %s, glom #%d, fly #%d', ...
%               ccoord(2),ccoord(1),r2_or_r4,glomnum_for_averagekernel, ...
%               compkernels(i).flynum));
%     end
%     
%     for i = 1:numel(compkernels)
%         set(0,'CurrentFigure',i);
%         ylim([ymin ymax]);
%         
%         if save_figures
%             savefig(sprintf('panoconv %d',i));
%         end
%     end
end

function tf = celleq(eqval,varargin)
    tf = cellfun(@(x)eq(x,eqval),varargin);
end

function cent = scalecent(cent,origimsize,newimsize)
%     fprintf('old cent: %f,%f\n',cent);
    cent = 1+(newimsize(:,[2 1])-1).*(cent-1)./(origimsize(:,[2 1])-1);
%     fprintf('new cent: %f,%f\n',cent);
end