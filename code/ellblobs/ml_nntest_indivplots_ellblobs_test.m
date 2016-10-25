function ml_nntest_indivplots_ellblobs_test
    doplot = true;
    dosave = true;
    ptrain = 0.4;
    sz = [16 10];
    
    suffix = '_offset';
    
    load(['ml_inputs_azimuths' suffix '.mat'],'t','r2_ind','r4_ind');
    
    load(['./ellblobs/ellblob' suffix '.mat']);
    
    ind = randperm(size(t,1));
    traini = ind(1:ptrain*length(ind));
    testi = ind(ptrain*length(ind)+1:end);
    
    if doplot
        figure(1);clf
        alsubplot(2,size(t,2),1,1);
    end
    
    mi = NaN(4,size(t,2));
    cols = 'brgm';
    mi(1,:)=dotest(x_im,t,traini,testi,1,doplot,cols);
    mi(2,:)=dotest(x_kern(:,r2_ind),t,traini,testi,2,doplot,cols);
    mi(3,:)=dotest(x_kern(:,r4_ind),t,traini,testi,3,doplot,cols);
    mi(4,:)=dotest(x_kern,t,traini,testi,4,doplot,cols);
    labels = {'views','R2','R4d','R2+R4d'};
    for i = 1:size(t,2)
        alsubplot(2,i)
        set(gca,'FontSize',8)
        colorbar(mi(:,i),cols,labels);
        ylim([0 1.6])
        if i==1
            ylabel('Mutual information (bits)')
        else
            set(gca,'YTick',[]);
        end
    end

    if dosave
        savefig('ell_colorbar',sz);
    end
end

function mi=dotest(x,t,traini,testi,figno,doplot,cols)
    nhidden = 10;
    nmibin = 100;
    iscirc = [true,false,false];
    actfun = 'linear';
    

    labels = {'Orientation','Size','Elevation'};
%     for i = 1:size(x,3)
    net = mlp(size(x,2),nhidden,size(t,2),actfun);
    options = zeros(1,18);
    % options(1) = 1;			% This provides display of error values.
    options(14) = 100;		% Number of training cycles. 
    net = netopt(net,options,x(traini,:),t(traini,:),'scg');
    y = mlpfwd(net,x(testi,:));

%     y = bsxfun(@minus,y_test,min(y_test));

%     if doplot
%         figure(1+(figno-1)*size(x,3));clf
%         alsubplot(1,size(y,2),1,1);
%     end
    mi = NaN(1,size(y,2));
    for j = 1:size(y,2)
        ut = unique(t(testi,j));
        [yy,yerr] = deal(NaN(size(ut)));
        mimat = NaN(size(y,1),length(ut));
        for k = 1:length(ut)
            cy = y(t(testi,j)==ut(k),j);
            mimat(1:length(cy),k) = cy;

%             if iscirc(j)
%                 yy(k) = mod((90/pi)*circ_mean((pi/90)*cy),180);
%                 yerr(k) = (90/pi)*circ_std((pi/90)*cy)./sqrt(length(cy));
%             else
                yy(k) = mean(cy);
                yerr(k) = std(cy)./sqrt(length(cy));
%             end
        end

        % calculate ideal observer metric
        [mi(j),p_smax]=ml_mi(mimat,nmibin);

        if doplot
            alsubplot(1,j)
            hold on
            set(gca,'FontSize',8);
            title(labels{j})
            if j==1
                ylabel('Network activation');
            end
            errorbar(ut,yy,yerr,cols(figno));
            xlim([0 1])
            ylim([0 1])
            if j>1
                set(gca,'YTick',[]);
            end
            
%             if j==2
%                 title(ttl)
%             end

%             lims = [min(ut) max(ut)];
%             xlim(lims)
% %                 ylim(lims)
% 
%             ylabel('Network activation');
%             xlabel(sprintf('%s (%.4f bits)',labels{j},mi))

%             alsubplot(2,j)
%             plot(ut,p_smax_test);

%             alsubplot(2,j)
%             errorbar(ut,yy_train,yerr_train,cols(j));
%             if j==2
%                 title('(ellipse test)')
%             end
%             lims = [min(ut) max(ut)];
%             xlim(lims)
% %                 ylim(lims)
% 
%             ylabel('Network activation');
%             xlabel(sprintf('%s (%.4f bits)',labels{j},mi_train))

%             alsubplot(4,j)
%             plot(ut,p_smax_train);
        end
    end
%     if doplot && dosave
%         savefig(sprintf('ml_ell2%s_%s',suffix,ttl));
%         close all
%     end
end