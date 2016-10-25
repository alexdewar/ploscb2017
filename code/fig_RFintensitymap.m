function fig_RFintensitymap
    r2_or_r4 = 'r4';

    load('vf_kernels.mat');
    main(eval(['vf_avkernels_' r2_or_r4]));
end

function main(kstruct)
    xlims = [-135 135];
    ylims = [-60 60];
    filtsize = 5;

    kstruct = kstruct(cell2mat({kstruct.isleft})); % & cell2mat({kstruct.glomnum})==12);
    kernels = {kstruct.k};
    imsz = size(kstruct(1).k);
    
%     indivgns = cell2mat({kstruct.glomnum});
%     gns = unique(indivgns);
%     shft = NaN(length(gns),2);
%     for i = 1:numel(gns)
%         avcent = mean(cell2mat({kstruct(indivgns==gns(i)).cent}'));
%         shft(i,:) = avcent; %-(imsz([2 1])/2);
%     end

    figure(2);clf
    hold on
    
    imsz = imsz-1;
    nkern = numel(kernels);
    alpha = 1/nkern;
    for i = 1:nkern
%         gni = gns==kstruct(i).glomnum;
%         curk = shiftmat(kernels{i},shft(gni,1),shft(gni,2));
        curk = kernels{i};

        [exx,exy] = bw2polygon(medfilt2(curk>0,filtsize*[1 1]));
        [exx,exy] = xytodeg(exx,exy);
        alfill(exx,exy,'r','linestyle','none','facealpha',alpha);
        line(exx,exy,'Color','r');

        [inx,iny] = bw2polygon(medfilt2(curk<0,filtsize*[1 1]));
        [inx,iny] = xytodeg(inx,iny);
        alfill(inx,iny,'b','linestyle','none','facealpha',alpha);
        line(inx,iny,'Color','b');
        
%         plot(shft(1),shft(2),'g+');
    end
    axis equal
    xlim(xlims);
    ylim(ylims);
%     keyboard

    function [x2,y2]=xytodeg(x,y)
        x2 = xlims(1) + range(xlims).*(x-1)/imsz(2);
        y2 = ylims(1) + range(ylims).*(imsz(1)-y+1)/imsz(1);
    end
end