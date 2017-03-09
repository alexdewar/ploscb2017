ymax = .3;

load([mfiledir '/ablate_data.mat']);

truesz = [36 28 14 42];

figure(1);clf
alsubplot(2,length(conds),1,1)
for i = 1:length(conds)
    alsubplot(1,i);
    plot(scores{i});
    hold on
    line(truesz(i)*[1 1],[0 ymax],'Color','b','LineStyle','--')
    
    ylabel('Total error')
    xlabel('Number of inputs')
    title(lconds{i});
    ylim([0 ymax]);
    xlim([1 length(scores{i})])
    
    alsubplot(2,i);
    bar(nanmean(inputscores{i}));
end