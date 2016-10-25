function panoconvs(dosave)
if ~nargin
    dosave = false;
end

dname = 'antoinestim';
dodiffs = false;

fulldname = fullfile(mfiledir,dname);
d = dir(fullfile(fulldname,'*.jpg'));

for i = 1:length(d)
    fname = fullfile(fulldname,d(i).name);
    fnnoext = d(i).name(1:end-4);

    figure(i);clf

    alsubplot(3 + 2*dodiffs,1,2,1);
    [vl,vr,ths,im]=vf_panoconv_average_centre_antoine(fname,'r2');
    meanactr2 = [mean(vl);mean(vr)]';
    plot(ths,meanactr2)
    title('R2');
    axis tight

    alsubplot(3,1);
    [vl,vr,ths]=vf_panoconv_average_centre_antoine(im,'r4');
    meanactr4 = [mean(vl);mean(vr)]';
    plot(ths,meanactr4)
    title('R4d')
    axis tight

    alsubplot(1,1)
    imagesc(im);
    colormap gray
    title(fnnoext,'Interpreter','none')
    axis off

    if dodiffs
        alsubplot(4,1)
        plot(ths(2:end),diff(meanactr2));
        title('R2 differential')
        axis tight

        alsubplot(5,1)
        plot(ths(2:end),diff(meanactr4));
        title('R4 differential')
        axis tight
    end
    
    if dosave
        savefig([dname '_' fnnoext],[10 10])
        close(i)
    end
end