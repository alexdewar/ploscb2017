%% TODO fix meddling and maybe add -ve vals for diffr2

function nsd_panoconvs_ridf_scores_r2(dosave)
    if ~nargin
        dosave = false;
    end
    
    doall = false;
    showprogbar = true;
    
    colcnt = [19 19];
    
    tcnt = [0,cumsum(colcnt)];
    ncol = length(colcnt);
    
    showrng = [45 135 225];
    im_size = [120 360];

    dname = '../../data/antoinestim';
    
    load('vf_kernels','vf_avkernels*');

    fulldname = fullfile(mfiledir,dname);
    if ~doall
        fulldname = [fulldname,'/touse'];
    end
    d = [dir(fullfile(fulldname,'*.jpg'));dir(fullfile(fulldname,'*.png'))];
    fname = sort({d.name});

    crng = round(showrng*im_size(2)/360);
    maxcol = max(colcnt);
    [diffr2,stdr2,eh_diffr2,retol] = deal(NaN(length(fname),1));
    eh_sig = zeros(size(diffr2));
    ims = ones(im_size(1)*size(diffr2,1),range(crng)+1,ncol);
    
    datafn = sprintf('%s/../dispatchrecap/panoconv.mat',mfiledir);
    if exist(datafn,'file')
        load(datafn);
    else
        if showprogbar
            startprogbar(1,length(fname));
        end
        
        ehf = fopen([mfiledir,'/dcp.txt'],'r');
        ehcell = textscan(ehf,'%d:%s\n');
        fclose(ehf);
        eh_fig = ehcell{1};
        eh_str = ehcell{2};

        for i = 1:length(fname)
            im = imread(fullfile(fulldname,fname{i}));
            if size (im,3)>1;
                im = rgb2gray(im);
            end
            im = im2double(im);

%             unitsperpix = 0.1*str2double(fname{i}(13))/str2double(fname{i}(9:11));
%             eh_diffr2(j,i) = unitsperpix*str2double(fname{i}(15:18));
%                 unitsperpix = eval(ehupp{ind});
%                 eh_diffr2(j,i) = unitsperpix*ehpx(ind);
            whcol = find(i <= tcnt,1)-1;
            ims((i-tcnt(whcol))*im_size(1)+(1:im_size(1)),:,whcol) = im(:,1+(crng(1):crng(3)));
            
            patt1 = im(:,crng(1):crng(2)-1);
            patt2 = im(:,crng(2):crng(3)-1);
%             figure(1);clf;
%             subplot(1,3,1);
%             imshow(patt1);
%             subplot(1,3,2);
%             imshow(patt2);
%             subplot(1,3,3);
%             imagesc(patt1 & patt2)
%             keyboard
            b_patt1 = ~patt1(:);
            b_patt2 = ~patt2(:);
            Q = sum(b_patt1 & b_patt2);
            R = sum(~b_patt1 & b_patt2);
%             T = sum(b_patt1 & ~b_patt2);
            retol(i) = Q/(Q+R);
%             fprintf('%s: %f, %f\n',fname{i},retol(i),Q/(Q+sum(~b_patt1 & b_patt2)))

            [diffr2(i),stdr2(i)] = vf_ridf(im,vf_avkernels_r2);

%                 jind = maxcol-j+1;
            eh_vals = eval(['[',eh_str{i},']']);
            eh_diffr2(i) = eh_vals(1);
            eh_sig(i) = eh_vals(2);
%             diffr2(i) = rr2(thr2==-90);
%             stdr2(i) = std(rr2)/sqrt(length(rr2));

%             sig(j,i) = str2double(fname{i}(6));

            if showprogbar && progbar
                return;
            end
        end
        save(datafn,'diffr2','stdr2','eh_diffr2','eh_fig','retol','ims','eh_sig');
    end % endhash
    
    
    %% plot figure
%     barsp = 0.1;
    
    sz = [19 10];
    
    xlo = 0; %-0.05;
    xhi = 0.4;
    xtick = 0.1;
    
    barsp = -xlo*im_size(1)/im_size(2);
    
%     imsc = 0.5;
%     
%     imsc = imsc*barsp;
    
%% proper fig code
    figure(1);clf
%     for i = 1:length(fname)
%         whcol = find(i <= tcnt,1)-1;
%         whrow = i-tcnt(whcol);
%         
%         alsubplot(ncol+1,2,[whcol whcol],1:2)
%         hold on
%         
% %         if whrow==1
% %             x = linspace(xlo,0,im_size(2));
% %             y = barsp*(0.5+linspace(maxcol,0,size(ims,1)));
% %             image(x,y,ims(:,:,whcol)*255);
% %             colormap gray
% %         end
%         
%         if eh_sig(i)==0
% %             barh(barsp*(maxcol-whrow+1),diffr2(i),0.8*barsp,'FaceColor','w');
%             bar(maxcol-whrow+1,diffr2(i),'FaceColor','w');
%         else
% %             barh(barsp*(maxcol-whrow+1),diffr2(i),0.8*barsp,'FaceColor',0.75*[1 1 1]);
%             bar(maxcol-whrow+1,diffr2(i),'FaceColor',0.75*[1 1 1]);
%         end
%         
%         errorbar(maxcol-whrow+1,diffr2(i),0,stdr2(i),'Color','k');
%         
% %         axis equal
% %         ylim(barsp*[0.25 1+maxcol])
%         ylim([xlo xhi])
% %         set(gca,'YTick',[]);
% %         set(gca,'XTick',0:xtick:xhi);
% %         xlabel('r.m.s. difference between 0^\circ and 90^\circ')
%     end
    
    %% temp fig code
%     figure(2);clf
%     bar([normalizevals(retol),normalizevals(abs(eh_diffr2)),normalizevals(abs(diffr2))])
% %     bar([retol,eh_diffr2,diffr2])
%     set(gca,'XTick',1:length(diffr2),'XTickLabel',eh_fig)
    
    %% stats
    sel = ~isnan(eh_diffr2);
    r2X = normalizevals(diffr2(sel));
    rolX = retol(sel);
    Y = abs(eh_diffr2(sel));
    [rhor2,pr2] = corr(r2X,Y,'type','Spearman');
    [rhorol,prol] = corr(rolX,Y,'type','Spearman');
    fprintf('R2: N = %d; rho = %f; p = %f\nROL: N = %d; rho = %f; p = %f\n', ...
            length(r2X),rhor2,pr2,length(r2X),rhorol,prol);
    
%     alsubplot(3,1)
    plot(r2X,Y,'kx')
%     hold on
%     r2c = r2Y-rhor2*X;
%     yl = get(gca,'YLim');
%     line(
    lsline
    xlabel('R2 difference')
    ylabel('Fly learning index')

%     alsubplot(3,2)
%     plot(rolX,Y,'kx')
%     lsline
%     xlabel('Retinal overlap')
    
        
    if dosave
        savefig('nsd_pattern_scatter',sz,'eps');
        
        fid = fopen(sprintf('%s/panoconv_%s.txt',mfiledir,mfilehash),'w');
        fprintf(fid,'Spearman''s correlation:\nrho = %f\np = %f\n',rhor2,pr2);
        fclose(fid);
    else
        dump2base(true)
    end
end

function [dr2,stdr2]=vf_ridf(im,kerns)
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
    
    acts0 = acts(:,1+size(rim,2)/2);
    rr2 = sqrt(mean(bsxfun(@minus,acts,acts0).^2));
    rr2 = [rr2,rr2(1)];
    
%     ths = linspace(-180,180,size(rim,2)+1);
    
    i90 = 1+size(rim,2)/4;
    dr2 = rr2(i90); % sign(mean(acts0)-mean(acts(:,i90)))*
    stdr2 = std(rr2)/sqrt(length(rr2));
end

% function horzerr(x,y,h,barsp)
%     for i = 1:numel(x)
%         line(x(i)+[0 h(i)],[y(i) y(i)],'Color','k')
%         line(x(i)+h(i)*[1 1],y(i)+barsp*[-0.2 0.2],'Color','k')
%     end
% end