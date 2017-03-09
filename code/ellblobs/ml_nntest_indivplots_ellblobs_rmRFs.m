function ml_nntest_indivplots_ellblobs_rmRFs
    doplot = false;
    dosave = false;
    
    suffix = '_offset';
    
    load(['ml_inputs_azimuths' suffix '.mat']);
    x_im_train = x_im(:,:,1);
    x_kern_train = x_kern(:,:,1);
%     t = t;
    
%     x_im = x_im(:,:,1);
%     x_kern = x_kern(:,:,1);
    
    load(['./ellblobs/ellblob' suffix '.mat']);
    
%     dotest(x_im_train,t,x_im,t,sprintf('Views (%dpx)',size(x_im,2)),1,doplot,dosave,suffix);
    dotest(x_kern_train(:,r2_ind),t,doplot,dosave,suffix);
%     dotest(x_kern_train(:,r4_ind),t,x_kern(:,r4_ind),t,'R4',3,doplot,dosave,suffix);
%     dotest(x_kern_train,t,x_kern,t,'R2+R4',4,doplot,dosave,suffix);
%     if dosave
%         system(sprintf('pdftk "./figures/ml_ell*(%04d).pdf" cat output ml_ell.pdf',fcnt));
%     end
end

function dotest(x,t,doplot,dosave,suffix)
    nhidden = 10;
    nmibin = 100;
    iscirc = [true,false,false];
    ptrain = 0.4;

    labels = {'Orient','Size','El'};
    cols = 'brg';
    options = zeros(1,18);
    % options(1) = 1;			% This provides display of error values.
    options(14) = 100;		% Number of training cycles. 
    cx = x;
    scores = NaN(size(x,2)/2+1,size(t,2));
    startprogbar(1,size(x,2)/2+1);
    for i = 0:size(x,2)/2
        if i~=0
            cind = [1:i-1, i+1:(2*i-1), (2*i+1):size(x,2)];
            cx = x(:,cind);
        end
        
        ind = randperm(size(t,1));
        traini = ind(1:ptrain*length(ind));
        testi = ind(ptrain*length(ind)+1:end);
        
        net = mlp(size(x,2)-2*(i~=0),nhidden,size(t,2),'linear');
        net = netopt(net,options,cx(traini,:),t(traini,:),'scg');
        y = mlpfwd(net,cx(testi,:));
        y = bsxfun(@minus,y,min(y));

        if doplot
            figure(i);clf
            alsubplot(1,size(y,2),1,1);
        end
        for j = 1:size(y,2)
            ut = unique(t(testi,j));
            [yy,yerr] = deal(NaN(size(ut)));
            mimat = NaN(size(y,1)./length(ut),length(ut));
            for k = 1:length(ut)
                cy = y(t(testi,j)==ut(k),j);
                mimat(1:length(cy),k) = cy;

                if iscirc(j)
                    yy(k) = mod((90/pi)*circ_mean((pi/90)*cy),180);
                    yerr(k) = (90/pi)*circ_std((pi/90)*cy)./sqrt(length(cy));
                else
                    yy(k) = mean(cy);
                    yerr(k) = std(cy)./sqrt(length(cy));
                end
            end

            % calculate ideal observer metric
            [scores(i+1,j),p_smax_test]=ml_mi(mimat,nmibin);
%             [mi_train,p_smax_train]=ml_mi(mimat_train,nmibin);

            if doplot
                alsubplot(1,j)
                errorbar(ut,yy,yerr,cols(j));
                if j==2
                    title(sprintf('missing %d',i))
                end

                lims = [min(ut) max(ut)];
                xlim(lims)
    %                 ylim(lims)

                ylabel('Network activation');
                xlabel(sprintf('%s (%.4f bits)',labels{j},scores(i,j)))
            end
        end
        
        if progbar
            return;
        end
    end
    if doplot && dosave
        savefig(sprintf('ml_ell%s_%s',suffix,ttl));
        close all
    end
    
    figure(1);clf
    bar(bsxfun(@minus,scores(2:end,:),scores(1,:)));
    keyboard
end