% UP BARSIMS COUNT

function ml_nntest_azimuthnet_multiaz(dosave,whnet,whtodo,dolabelerrs,whcvar,figno,fixedv,fixedval)
%     if nargin < 2
%         whnet = 1;
%     end
%     if nargin < 1
%         dosave = false;
%     end
    
    netsuffixes = { 'elaz', 'orsiel+az', 'orsiel', 'elaz+or' };
    switch whnet
        case 1
            nvar = 2;
            vnames = {'el','az'};
        case 4
            nvar = 2;
            vnames = {'el','az','or'};
        otherwise
            nvar = 3;
            vnames = {'or','si','el','az'};
    end
    if nargin < 4
        dolabelerrs = false;
    end
    if nargin < 3
        whtodo = 1:nvar;
    else
        if ischar(whtodo)
            close all
            switch whtodo
                case 'perm'
                    for i = 1:nvar
                        ml_nntest_azimuthnet_multiaz(dosave,whnet,i,dolabelerrs,[],i);
                    end
                    ml_nntest_azimuthnet_multiaz(dosave,whnet,1:nvar,dolabelerrs,[],nvar+1);
                    return
                case 'all'
                    whtodo = 1:nvar;
            end
        end
    end
    
    if nargin < 8
        fixedv = [];
    end
    dofixedv = ~isempty(fixedv);
    
    if nargin < 6 || isempty(figno)
        figno = 1;
    end
    if nargin < 5
        whcvar = [];
    elseif ~isempty(whcvar)
        whtodo = whtodo(whtodo~=whcvar);
    end
    
    doonelayer = false;
    doloadpreprocess = false;
    doplot = true;
    
    figsz = [20 10];
    docombine = isempty(whcvar);
    dname = fullfile(mfiledir,'../drosodata/ANNs');
    prefix = ['ellblob_' netsuffixes{whnet}];
    ptrain = 0.4;

    r2_ind = 1:28;
    r4_ind = 28+(1:14);
    
    if whnet==2
        nvar = 4;
    end
    tdtf = false(1,nvar);
    tdtf(whtodo) = true;
    tdtfstr = sprintf('%d',tdtf);
    
    fprintf('========= doing %s (%s) =========\n',netsuffixes{whnet},tdtfstr)
    
%     if whnet==4
%         whnet = 1;
%     end
%     if whnet~=2
%         docombine = false;
%     end
    if docombine
        paramstr = '';
    else
        paramstr = ['_by_' vnames{whcvar}];
    end
    if dofixedv
        if ischar(fixedval)
            paramstr = sprintf('%s_%s=S[',paramstr,vnames{fixedv});
            paramstr = [paramstr, urlencode(fixedval), ']'];
        else
            paramstr = sprintf('%s_%s=%.2f',paramstr,vnames{fixedv},fixedval);
        end
    end
    if doonelayer
        paramstr = [paramstr '_1LAYER'];
        warning('DOING ONE LAYER NET')
    end
    
    d = dir(sprintf('%s/%s*.mat',dname,prefix));
    if length(d)~=1
        error('found %d ellblob files, expected 1',length(d))
    end
    fname_ellblob = d(1).name;
    fname_preprocess = fullfile(dname,sprintf('fig_%s_td%s_%s%s.mat',fname_ellblob,tdtfstr,prefix,paramstr));
    if doloadpreprocess && exist(fname_preprocess,'file')
        load(fname_preprocess);
    else
        load(fullfile(dname,fname_ellblob));
        
        switch whnet
            case 1
                t = [el, az];
            case 4
                t = [el, az, orients];
            otherwise
                t = [orients, areas/1000, el, az];
        end
        
        if dofixedv
            if ischar(fixedval)
                ut = unique(t(:,fixedv));
                fixedval = eval(['ut(' fixedval ')']);
            end
            
%             usz = unique(t(:,2));
            
            fvsel = t(:,fixedv)==fixedval; %& t(:,2)==usz(11);
            t = t(fvsel,:);
            x_im = x_im(fvsel,:);
            x_kern = x_kern(fvsel,:);
            
            
            fprintf('fixing %s at %.2f\n',vnames{fixedv},fixedval)
        end
        
        if docombine
            ncvar = 1;
            ucvar = [];
        else
            cvar = t(:,whcvar);
            ucvar = unique(cvar);
            ncvar = length(ucvar);
        end

        t = t(:,whtodo);
        
        ndata = size(t,1)/ncvar;
        
        ind = randperm(ndata);
        ntrain = round(ptrain*ndata);
        traini = ind(1:ntrain);
        testi = ind(ntrain+1:end);
        
        [meanval,errval] = deal(NaN(size(t,2),4,ncvar));
        ndataeach = length(unique(t(:,1)));
        [ut,yy,yerr] = deal(NaN(size(t,2),ndataeach,4,ncvar));
        
        if whnet==1 || whnet==4
            [ut(:,:,1),yy(:,:,1),yerr(:,:,1),meanval(:,1),errval(:,1)]=dotest(x_im,t,traini,testi,sprintf('Views (%dpx)',size(x_im,2)),ndataeach,doonelayer);
            [ut(:,:,2),yy(:,:,2),yerr(:,:,2),meanval(:,2),errval(:,2)]=dotest(x_kern(:,r2_ind),t,traini,testi,'R2',ndataeach,doonelayer);
            [ut(:,:,3),yy(:,:,3),yerr(:,:,3),meanval(:,3),errval(:,3)]=dotest(x_kern(:,r4_ind),t,traini,testi,'R4',ndataeach,doonelayer);
            [ut(:,:,4),yy(:,:,4),yerr(:,:,4),meanval(:,4),errval(:,4)]=dotest(x_kern,t,traini,testi,'R2+R4',ndataeach,doonelayer);
        else
            if whnet==2
                fprintf('===== across %s =====\n',joinstr(', ',vnames{~ismember(1:4,[whtodo,whcvar,fixedv])}))
            end
            if ~docombine
                fprintf('===== grouped by %s =====\n',vnames{whcvar})
            end
            for i = 1:ncvar
                if docombine
                    if whnet~=2
                        fprintf('===== all vars combined =====\n')
                    end
                    sel = true(size(t,1),1);
                else
                    fprintf('===== %s %d/%d =====\n',vnames{whcvar},i,ncvar)
                    sel = cvar==ucvar(i);
                end
                [ut(:,:,1,i),yy(:,:,1,i),yerr(:,:,1,i),meanval(:,1,i),errval(:,1,i)]=dotest(x_im(sel,:),t(sel,:),traini,testi,sprintf('Views (%dpx)',size(x_im,2)),ndataeach,doonelayer);
                [ut(:,:,2,i),yy(:,:,2,i),yerr(:,:,2,i),meanval(:,2,i),errval(:,2,i)]=dotest(x_kern(sel,r2_ind),t(sel,:),traini,testi,'R2',ndataeach,doonelayer);
                [ut(:,:,3,i),yy(:,:,3,i),yerr(:,:,3,i),meanval(:,3,i),errval(:,3,i)]=dotest(x_kern(sel,r4_ind),t(sel,:),traini,testi,'R4',ndataeach,doonelayer);
                [ut(:,:,4,i),yy(:,:,4,i),yerr(:,:,4,i),meanval(:,4,i),errval(:,4,i)]=dotest(x_kern(sel,:),t(sel,:),traini,testi,'R2+R4',ndataeach,doonelayer);
                disp('.')
            end
        end
        
        if doloadpreprocess
            savemeta(fname_preprocess,'ut','yy','yerr','meanval','errval','fname_ellblob','t','ucvar','ncvar','whcvar','pgeb')
        end
    end
    
    if doplot
        figure(figno);clf
%         set(gca,'FontSize',8);
%         alsubplot(size(t,2),4,1,1);
        alsubplot(max(ncvar,2),size(t,2),1,1);
        
        if ncvar>1
            figure(100+figno);clf
            alsubplot(ncvar,size(t,2),1,1);
        end
    end
    
    cols = 'brgm';
    labels = {'views','R2','R4d','R2+R4d'};
    
    if whnet==1 || whnet==4
        tlabels = {'Elevation (deg)','Azimuth (deg)'};
    else
        tlabels = {'Orientation (deg)','Size (x10^3 deg^2)','Elevation (deg)','Azimuth (deg)'};
    end
    
    for i = 1:ncvar
        for j = 1:size(t,2)
            if ncvar==1
                alsubplot(1,j)
            else
                figure(figno)
                alsubplot(i,j)
            end
            set(gca,'FontSize',8)
            hold on

            for k = 1:4
                plot(ut(j,:,k,i),yy(j,:,k,i),cols(k))
                ciplot(yy(j,:,k,i)-yerr(j,:,k,i),yy(j,:,k,i)+yerr(j,:,k,i),ut(j,:,k,i),cols(k),'EdgeColor','none')
            end
            lims = [min(min(ut(j,:,:,i),[],3)) max(max(ut(j,:,:,i),[],3))];
            xlim(lims)
            ylim(lims)
            axis square

%             switch j
%                 case 1
%                     if naz==1
%                         ylabel('Network activation');
%                     end
%                     set(gca,'XTick',0:10:90,'YTick',0:10:90);
%                 case 3
%                     set(gca,'XTick',-60:20:60,'YTick',-60:20:60);
%             end

            if i==ncvar
                xlabel(tlabels{whtodo(j)})
            else
                set(gca,'XTick',[])
            end
            if whnet~=1 && whnet~=4 && ~docombine && j==ceil(size(t,2)/2)
                title(sprintf('%s=%.2f',vnames{whcvar},ucvar(i)))
            end
            if j==1
                if ncvar==1
                    ylabel('Mean square error')
                end
%             else
%                 set(gca,'YTick',[]);
            end
            

    %             fprintf('%s (%.4f bits)\n',labels{j},mi(j))

    %             alsubplot(2,j)
    %             hold on
    %             bar(figno,mi,cols(j))

            if ncvar==1
                alsubplot(2,j)
            else
                figure(100+figno)
                alsubplot(i,j)
                
                if ~docombine && j==ceil(size(t,2)/2)
                    title(sprintf('%s=%.2f',vnames{whcvar},ucvar(i)))
                end
            end
            set(gca,'FontSize',8)
            hold on
            if i==ncvar
                alcolorbar(meanval(j,:,i),cols,labels);
            else
                alcolorbar(meanval(j,:,i),cols,[]);
            end
            if ncvar>1 || dolabelerrs
                for k = 1:size(meanval,2)
                    text(k,meanval(j,k,i)+errval(j,k,i),sprintf(' %.3f',meanval(j,k,i)),'Rot',90);
                end
            end
            errorbar(1:4,meanval(j,:,i),zeros(1,4),errval(j,:,i),'k','LineStyle','none');
    %         alpha(.5)
            ylim([0 .1])
            axis square

    %         alsubplot(3,i)
    %         set(gca,'FontSize',8)
    %         colorbar(corrval(:,i,2),cols,labels);
    % %         alpha(.5)
    %         ylim([0 .2])
    %         if i==1
    %             ylabel('Mean square error')
    %         else
    %             set(gca,'YTick',[]);
    %         end
    %         axis square
        end
    end
%     end

%     dump2base(true)
    
    if doplot && dosave
        fstr = sprintf('ml_%s_%s%s',netsuffixes{whnet},tdtfstr,paramstr);
        if ncvar==1
            savefig(fstr,figsz,'pdf');
        else
            figsz(2) = figsz(2)*ncvar/2;
            savefig([fstr '_mse'],figsz,'pdf');
            figure(figno)
            savefig([fstr '_plot'],figsz,'pdf');
        end
    end
end

function [ut,yy,yerr,meanval,errval]=dotest(x,t,traini,testi,ttl,ndataeach,doonelayer)
    disp(ttl)
    
    rngt = range(t);
    
    mins = min(t);
%     maxs = max(t);
    normt = bsxfun(@rdivide,bsxfun(@minus,t,mins),rngt);

    nhidden = 10;
    nbarsims = 1;

    if doonelayer
        nbarsims = 1;
    end
    errors = NaN(nbarsims,size(t,2));
    for i = 1:nbarsims
        if doonelayer
            xbias = [ones(size(x,1),1), x]';
            w = normt'*pinv(xbias);
            y = (w*xbias)';
            errors(i,:) = mean((normt-y).^2);
        else
            net = mlp(size(x,2),nhidden,size(t,2),'linear');
            options = zeros(1,18);
            % options(1) = 1;			% This provides display of error values.
            options(14) = 100;		% Number of training cycles. 
            net = netopt(net,options,x(traini,:),normt(traini,:),'scg');
            y = mlpfwd(net,x(testi,:));
            errors(i,:) = mean((normt(testi,:)-y).^2);
        end
        
        y = bsxfun(@plus,mins,bsxfun(@times,rngt,y));
    end
    
    [ut,yy,yerr] = deal(NaN(size(t,2),ndataeach));
    for i = 1:size(y,2)
        ut(i,:) = unique(t(testi,i));
%         [yy,yerr] = deal(NaN(size(ut)));
        for j = 1:length(ut)
            cy = y(t(testi,i)==ut(i,j),i);
            
            if any(isnan(cy))
                error('nans in cy!')
            end

%             if iscirc(j)
%                 yy(k) = mod((90/pi)*circ_mean((pi/90)*cy),180);
%                 yerr(k) = (90/pi)*circ_std((pi/90)*cy)./sqrt(length(cy));
%             else
                yy(i,j) = mean(cy);
                yerr(i,j) = std(cy)./sqrt(length(cy));
%             end
        end

%         errval(1,i,1) = mean((t(testi,i)-y(:,i)).^2);
%         sel = t(testi) >= 0.25 & t(testi) <= 0.75;
%         errval(1,i,2) = mean((t(testi(sel),i)-y(sel,i)).^2);

%         if doplot
% 
%         end
    end
    disp('.');
    
    meanval = mean(errors,1)';
    errval = stderr(errors,1)';
end