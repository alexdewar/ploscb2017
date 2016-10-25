function panoconvs_ridf_scores_r2_debug(dosave)
    if ~nargin
        dosave = false;
    end
    
    doall = false;
    showprogbar = true;
    
    ncol = 2;
    colcnt = [19 13];
    
    showrng = [45 225];
    im_size = [120 360];

    dname = 'antoinestim';
    
    load('vf_kernels','vf_avkernels*');

    fulldname = fullfile(mfiledir,dname);
    if ~doall
        fulldname = [fulldname,'/touse'];
    end
    d = [dir(fullfile(fulldname,'*.jpg'));dir(fullfile(fulldname,'*.png'))];
    fname = sort({d.name});

    crng = round(showrng*im_size(2)/360);
    % ceil(length(fname)/ncol)
    maxcol = max(colcnt);
    [diffr2,stdr2,eh_diffr2] = deal(NaN(maxcol,ncol));
    sig = zeros(size(diffr2));
    ims = ones(im_size(1)*size(diffr2,1),range(crng)+1,ncol);
    
    datafn = sprintf('%s/panoconv_%s.mat',mfiledir,mfilehash);
    if exist(datafn,'file')
        load(datafn);
    else
        if showprogbar
            startprogbar(1,length(fname));
        end

        toti = 0;
        for i = 1:ncol
            for j = 1:colcnt(i)
                ind = toti+j;
                if ind>numel(fname)
                    break;
                end

                im = imread(fullfile(fulldname,fname{ind}));
                if size (im,3)>1;
                    im = rgb2gray(im);
                end
                im = im2double(im);
                
                unitsperpix = 0.1*str2double(fname{ind}(13))/str2double(fname{ind}(9:11));
                eh_diffr2(j,i) = unitsperpix*str2double(fname{ind}(15:18));
                
                ims((j-1)*im_size(1)+(1:im_size(1)),:,i) = im(:,1+(crng(1):crng(2)));

                [rr2,thr2]=vf_ridf(im,vf_avkernels_r2);
    
%                 jind = maxcol-j+1;

                diffr2(j,i) = rr2(thr2==-90);
                stdr2(j,i) = std(rr2)/sqrt(length(rr2));

                sig(j,i) = str2double(fname{ind}(6));

                if showprogbar && progbar
                    return;
                end
            end
            toti = toti+colcnt(i);
        end
        save(datafn,'diffr2','stdr2','eh_diffr2','sig','ims');
    end % endhash
    
    
    %% plot figure
%     barsp = 0.1;
    
    sz = [19 10];
    
    xlo = -0.05;
    xhi = 0.4;
    xtick = 0.1;
    
    barsp = -xlo*im_size(1)/im_size(2);
    
    for figno = 1:2
        if figno==1
            data = diffr2;
            dataerr = stdr2;
        else
            data = eh_diffr2;
            dataerr = zeros(size(eh_diffr2));
        end
        
        figure(figno);clf

        for i = 1:ncol
            subplot(1,ncol,i)
            hold on

            x = linspace(xlo,0,im_size(2));
            y = barsp*(0.5+linspace(maxcol,0,size(ims,1)));
            image(x,y,ims(:,:,i)*255);
            colormap gray

            for j = 1:colcnt(i)
                if sig(j,i)==0
                    barh(barsp*(maxcol-j+1),data(j,i),0.8*barsp,'FaceColor','w');
                else
                    barh(barsp*(maxcol-j+1),data(j,i),0.8*barsp,'FaceColor',0.75*[1 1 1]);
                end
            end

            horzerr(data(:,i),barsp*(size(data,1):-1:1),dataerr(:,i),barsp);

            axis equal
            ylim(barsp*[0.25 1+maxcol])
            xlim([xlo xhi])
            set(gca,'YTick',[]);
            set(gca,'XTick',0:xtick:xhi);
    %         xlabel('r.m.s. difference between 0^\circ and 90^\circ')
        end
    end
    
    X = diffr2(:);
    Y = eh_diffr2(:);
    sel = ~isnan(Y);
    X = X(sel);
    Y = Y(sel);
    [rho,p] = corr(X,Y,'type','Spearman');
    fprintf('Spearson''s correlation:\nrho = %f\np = %f\n',rho,p);
    
    if dosave
        savefig('pattern_ridf',sz);
        
        fid = fopen(sprintf('%s/panoconv_%s.txt',mfiledir,mfilehash),'w');
        fprintf(fid,'Spearson''s correlation:\nrho = %f\np = %f\n',rho,p);
        fclose(fid);
    else
        dump2base(true)
    end
    
end

function [ridf,ths]=vf_ridf(im,kerns)
    kerns = cell2mat(shiftdim({kerns.k},-1));
    
    ksz = size(kerns);
    rim = imresize(im,[ksz(1),ceil(360*ksz(2)/270)]);
    xoff = (size(rim,2)-ksz(2))/2;
    acts = NaN(ksz(3),size(rim,2));
    for i = 1:size(rim,2)
        crim = circshift(rim,[0 i-1-ceil(size(rim,2)/2)]);
        
        acts(:,i) = shiftdim(sum(sum(bsxfun(@times,crim(:,xoff+(1:ksz(2))),kerns)),2));
    end
    
%     for i = 1:size(im,2)
%         [acts(:,i),kerns] = getneuronactivations(circshift(im,[0 i-1-ceil(size(im,2)/2)]),kerns);
%     end
    
    ridf = sqrt(mean(bsxfun(@minus,acts,acts(:,1+size(rim,2)/2)).^2));
    ridf = [ridf,ridf(1)];
    
    ths = linspace(-180,180,size(rim,2)+1);
end

function horzerr(x,y,h,barsp)
    for i = 1:numel(x)
        line(x(i)+[0 h(i)],[y(i) y(i)],'Color','k')
        line(x(i)+h(i)*[1 1],y(i)+barsp*[-0.2 0.2],'Color','k')
    end
end