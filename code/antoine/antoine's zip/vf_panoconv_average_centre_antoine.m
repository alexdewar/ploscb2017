function [vals_l,vals_r,ths,im]=vf_panoconv_average_centre_antoine(image_filename, r2_or_r4)
    load('vf_kernels.mat');
        avkernels = eval(['vf_avkernels_' r2_or_r4]);
    for i=1: length(avkernels)/2; %For each kernels
     glomnum_for_averagekernel = i; 
    [val_l,val_r,ths,im]=mainfun(image_filename,thresh,vf_kernels_r2,vf_kernels_r4,vf_avkernels_r2,vf_avkernels_r4,r2_or_r4,glomnum_for_averagekernel);
    vals_l (i,:) = val_l;
    vals_r (i,:) = val_r;
    end
end

function [vals_l,vals_r,ths,im]=mainfun(im,thresh,vf_kernels_r2,vf_kernels_r4,vf_avkernels_r2,vf_avkernels_r4,r2_or_r4,glomnum_for_averagekernel)
%     r2_or_r4 = 'r2';
    

%     use_left_version_of_kernel = true;
%     image_filename = 'steps.jpg';
    save_figures = false;
    show_figures = false;

    thresh = 0.25;
    krng = [-60 60; -135 135]';
    irng = [-60 60; -180 180]';
    
    kernels = eval(['vf_kernels_' r2_or_r4]);
    avkernels = eval(['vf_avkernels_' r2_or_r4]);

    if ischar(im)
        im = imread(im);
    end
    if size (im,3)>1;
        im = rgb2gray(im);
    end
    im = im2double(im);
    
    ifov = range(irng);
    kfov = [ifov(1),range(krng(:,2))]; % squish kernel vertically
    
% for i=1: length(avkernels)/2; %For each kernels
%     glomnum_for_averagekernel = 2;%i; 
    
    avkernels = avkernels(celleq(glomnum_for_averagekernel,avkernels.glomnum));
    % kernels to use to calculate where centre should be for average kernel
    compkernels = kernels(celleq(glomnum_for_averagekernel,kernels.glomnum));
    
    
    isleft = cell2mat({avkernels.isleft});
%     val_l =[]; val_r=[];
    vals_l =dokernel(avkernels(isleft),true);
    vals_r =dokernel(avkernels(~isleft),false);
    
%     vals_l (i,:) = val_l;
%     vals_r (i,:) = val_r;
%  end

    function vals=dokernel(kernel,isleft)
%         figure(10*isleft+1);clf
%         showkernel(kernel);
        
        %% resize kernel
        kern = (sign(kernel.k)+1)/2; % set kernel to be between 0 and 1

        % calculate new size
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
        avkcent = scalecent(kernel.cent,origksize,size(kern));
%         figure(10*isleft+2);clf
%         showkernel(kern,avkcent);
        
        kcents = cell2mat({compkernels(celleq(isleft,compkernels.isleft)).cent}');
        kcents = scalecent(kcents,repmat(origksize,size(kcents,1)), ...
                                  repmat(size(kern),size(kcents,1)));
        averageofcents = mean(kcents);

%         return
        ths = linspace(irng(1,2),irng(2,2),length(vals));
        ymax = 0;
        ymin = inf;
%         if isleft
%             startprogbar(20,2*length(vals),'Doing convolutions...');
%         end
    %     for i = 1:numel(compkernels)
            coff = round(averageofcents-avkcent);
            if coff(2) < 0
                ckern = [kern(1-coff(2):end,:);zeros(-coff(2),size(kern,2))];
            else
                ckern = [zeros(coff(2),size(kern,2));kern(1:end-coff(2),:)];
            end
            if coff(1) < 0
                ckern = [ckern(:,1-coff(1):end) zeros(size(ckern,1),-coff(1))];
            else
                ckern = [zeros(size(ckern,1),coff(1)) ckern(:,1:end-coff(1))];
            end

    %         figure(1);clf
    %         subplot(1,3,1);
    %         showkernel(imresize(compkernels(i).k,szfac),kcents(i,:));
    %         subplot(1,3,2);
    %         showkernel(kern,avkcent);
    %         subplot(1,3,3);
    %         showkernel(ckern,kcents(i,:));
    %         keyboard

            for j = 1:length(vals)
                vals(j) = sum(sum(ckern.*im2(:,j-1+(1:size(ckern,2)))));

%                 if progbar
%                     return;
%                 end
            end

            if show_figures
                maxval = max(vals);
                if maxval > ymax
                    ymax = maxval;
                end
                minval = min(vals);
                if minval < ymin
                    ymin = minval;
                end


                figure(isleft+1);clf
                plot(ths,vals);
                axis tight

                ksz = size(kern);
                prop = (averageofcents([2 1])-1)./ksz;
                ccoord = prop.*ifov+irng(1,:);
                if isleft
                    lrstr = 'left RF';
                else
                    lrstr = 'right RF';
                end
                title(sprintf('Centre at (%.2fdeg, %.2fdeg); %s, glom #%d, %s', ...
                      ccoord(2),ccoord(1),r2_or_r4,glomnum_for_averagekernel,lrstr));
        %     end
    %             if isleft
    %                 return
    %             end
        %     for i = 1:numel(compkernels)
                set(0,'CurrentFigure',1+isleft);
                ylim([ymin ymax]);

                if save_figures
                    savefig([im(1:end-4) ' ' lrstr]);
                end
            end
    %     end
    end
end

function tf = celleq(eqval,varargin)
    tf = cellfun(@(x)eq(x,eqval),varargin);
end

function cent = scalecent(cent,origimsize,newimsize)
%     fprintf('old cent: %f,%f\n',cent);
    cent = 1+(newimsize(:,[2 1])-1).*(cent-1)./(origimsize(:,[2 1])-1);
%     fprintf('new cent: %f,%f\n',cent);
end