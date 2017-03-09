function fig_RFintensitymap_3d
    r2_or_r4 = 'r2';

    load('vf_kernels.mat');
    main(eval(['vf_avkernels_' r2_or_r4]));
end

function main(kstruct)
    maxth = -inf;
    minth = inf;

    xlims = [-135 135];
    ylims = [-60 60];
    filtsize = 5;
    
    kstruct = kstruct(cell2mat({kstruct.isleft}));% & cell2mat({kstruct.glomnum})==8);
    kernels = {kstruct.k};
    imsz = size(kstruct(1).k);
    
    cylr = 1;
    xlimsr = pi*xlims/180;
    xrng = range(xlimsr);
    yrng = range(pi*ylims/180);
    cylht = yrng;
    
    figure(1);clf
    hold on
    
    nkern = numel(kernels);
    for i = 1:nkern
        curk = kernels{i};
        draw3d(curk>0,'r');
        draw3d(curk<0,'b');
    end

    lwid = 1;
    
    ang = pi*135/180;
    ths = pi+linspace(-ang,ang,1000);
    [linex,liney] = pol2cart(ths,cylr);
    linez = zeros(size(linex));

    line(linex,liney,linez,'Color','k','LineWidth',lwid);
    line(linex,liney,linez+cylht,'Color','k','LineWidth',lwid);
    line([0 0],-cylr*[1 1],[0 cylht],'Color','k','LineWidth',lwid);
    line([0 0],cylr*[1 1],[0 cylht],'Color','k','LineWidth',lwid);
    
    axis equal off
    xlim(2*cylr*[-1 1]);
    ylim(2*cylr*[-1 1]);
    
    set(gca,'CameraPosition',[27.3036 0 15.7676],'CameraViewAngle',7.55833);
    
%     [cx,cy] = pol2cart(pi+[-135 135]*pi/180,cylr);
%     plot(cx,cy,'g+','markersize',5)
    
%     savefig('drum',[15 15],'png');
    
%     fprintf('min: %f\nmax: %f\n',minth,maxth);

    plot3(0,0,cylht/2,'g+')

%     keyboard
    
    function draw3d(cim,col)
        [x,y] = bw2polygon(medfilt2(cim,filtsize*[1 1]));
        
        th = xrng*(imsz(2)-x+1)/(imsz(2)-1);
        [X3,Y3] = pol2cart(th,cylr);
        Z3 = yrng*(imsz(1)-y)/(imsz(1)-1);
        
%         cminth = min(th);
%         if cminth < minth
%             minth = cminth;
%         end
%         cmaxth = max(th);
%         if cmaxth > maxth
%             maxth = cmaxth;
%         end
%         
%         fprintf('ang range: %f\n',circ_rad2ang(range(th)));
        nans = [0;find(isnan(x))];
        for n = 2:length(nans)
            ind = nans(n-1)+1:nans(n)-1;
%             fprintf('%d gap: x=%f; y=%f; z=%f\n',n-1,max(circdiff(X3(ind))),max(circdiff(Y3(ind))),max(circdiff(Z3(ind))));
            
            h=fill3(X3(ind),Y3(ind),Z3(ind),col,'edgecolor',col,'facealpha',1/nkern);
%             keyboard
        end
    end
end