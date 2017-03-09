function ml_nntest_indivplots_ellblobs_trainellipse
    doplot = true;
    dosave = false;
    
    suffix = '_offset';
    
    load(['ml_inputs_azimuths' suffix '.mat']);
    x_im_train = x_im(:,:,1);
    x_kern_train = x_kern(:,:,1);
%     traint = t;
    
    load(['./ellblobs/ellblob' suffix '.mat']);
%     x_im_train = x_im(:,:,1);
%     x_kern_test = x_kern(:,:,1);
    
    dotest(x_im_train,t,x_im,t,sprintf('Views (%dpx)',size(x_im,2)),1,doplot,dosave,suffix);
    dotest(x_kern_train(:,r2_ind),t,x_kern(:,r2_ind),t,'R2',2,doplot,dosave,suffix);
    dotest(x_kern_train(:,r4_ind),t,x_kern(:,r4_ind),t,'R4',3,doplot,dosave,suffix);
    dotest(x_kern_train,t,x_kern,t,'R2+R4',4,doplot,dosave,suffix);
%     if dosave
%         system(sprintf('pdftk "./figures/ml_ell*(%04d).pdf" cat output ml_ell.pdf',fcnt));
%     end
end

function dotest(trainx,traint,testx,testt,ttl,figno,doplot,dosave,suffix)
    nhidden = 10;
    nmibin = 100;
    iscirc = [true,false,false];
    ptrain = 0.4;
    actfun = 'linear';
    
    ind = randperm(size(traint,1));
%     traini = ind;
%     testi = ind;
    traini = ind(1:ptrain*length(ind));
    testi = ind(ptrain*length(ind)+1:end);

    labels = {'Orient','Size','El'};
    cols = 'brg';
%     for i = 1:size(trainx,3)
    net_test = mlp(size(trainx,2),nhidden,size(traint,2),actfun);
    options = zeros(1,18);
    % options(1) = 1;			% This provides display of error values.
    options(14) = 100;		% Number of training cycles. 
    net_test = netopt(net_test,options,testx(traini,:),traint(traini,:),'scg');
    y_test = mlpfwd(net_test,testx(testi,:));

    net_train = mlp(size(trainx,2),nhidden,size(traint,2),actfun);
    options = zeros(1,18);
    % options(1) = 1;			% This provides display of error values.
    options(14) = 100;		% Number of training cycles. 
    net_train = netopt(net_train,options,trainx(traini,:),traint(traini,:),'scg');
    y_train = mlpfwd(net_train,trainx(testi,:));
    y_train = bsxfun(@minus,y_train,min(y_train));

    y_test = bsxfun(@minus,y_test,min(y_test));

    if doplot
        figure(1+(figno-1)*size(trainx,3));clf
        alsubplot(2,size(y_test,2),1,1);
    end
    for j = 1:size(y_test,2)
        ut = unique(traint(testi,j));
        [yy_test,yerr_test,yy_train,yerr_train] = deal(NaN(size(ut)));
        [mimat_test,mimat_train] = deal(NaN(size(y_test,1),length(ut)));
        for k = 1:length(ut)
            cy_test = y_test(traint(testi,j)==ut(k),j);
            mimat_test(1:length(cy_test),k) = cy_test;
            cy_train = y_train(traint(testi,j)==ut(k),j);
            mimat_train(1:length(cy_train),k) = cy_train;

%             if iscirc(j)
%                 yy_test(k) = mod((90/pi)*circ_mean((pi/90)*cy_test),180);
%                 yerr_test(k) = (90/pi)*circ_std((pi/90)*cy_test)./sqrt(length(cy_test));
% 
%                 yy_train(k) = mod((90/pi)*circ_mean((pi/90)*cy_train),180);
%                 yerr_train(k) = (90/pi)*circ_std((pi/90)*cy_train)./sqrt(length(cy_train));
%             else
                yy_test(k) = mean(cy_test);
                yerr_test(k) = std(cy_test)./sqrt(length(cy_test));

                yy_train(k) = mean(cy_train);
                yerr_train(k) = std(cy_train)./sqrt(length(cy_train));
%             end
        end

        % calculate ideal observer metric
        [mi_test,p_smax_test]=ml_mi(mimat_test,nmibin);
        [mi_train,p_smax_train]=ml_mi(mimat_train,nmibin);

        if doplot
            alsubplot(1,j)
            errorbar(ut,yy_test,yerr_test,cols(j));
            if j==2
                title(ttl)
            end

            lims = [min(ut) max(ut)];
            xlim(lims)
%             ylim(lims)

            ylabel('Network activation');
            xlabel(sprintf('%s (%.4f bits)',labels{j},mi_test))

%             alsubplot(2,j)
%             plot(ut,p_smax_test);

            alsubplot(2,j)
            errorbar(ut,yy_train,yerr_train,cols(j));
            if j==2
                title('(ellipse test)')
            end
            lims = [min(ut) max(ut)];
            xlim(lims)
%             ylim(lims)

            ylabel('Network activation');
            xlabel(sprintf('%s (%.4f bits)',labels{j},mi_train))

%             alsubplot(4,j)
%             plot(ut,p_smax_train);
        end
    end
    if doplot && dosave
        savefig(sprintf('ml_ell%s_%s',suffix,ttl));
        close all
    end
end