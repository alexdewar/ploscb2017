clear

dname = fullfile(mfiledir,'../../../data/rx_neurons/oldpaths20150219');

d = dir(fullfile(dname,'*.mat'));

errcnt = 0;
queues = {};
viewtypes = {'hires','lores','R2','R4','Rx'};
lenoff = 20;

[rts,time,whq,whvt] = deal(NaN(length(d),1));
for i = 1:length(d)
    fname = fullfile(dname,d(i).name);
    if ~all(varsinmatfile(fname,'runtime','metadata'))
        errcnt = errcnt+1;
        continue;
    end
    load(fname,'runtime','metadata');
    
    loff=length(fname)-lenoff;
    vt = fname(find(fname(1:loff)=='_',1,'last')+1:loff);
    whvt(i) = find(strcmpi(viewtypes,vt));
    
    rts(i) = runtime;
    time(i) = datenum(metadata.datetime)-runtime;
    
    curq = find(strcmp(metadata.sge_queue,queues));
    if isempty(curq)
        queues{end+1} = metadata.sge_queue;
        whq(i) = length(queues);
    else
        whq(i) = curq;
    end
end
if errcnt>0
    error('errors were found!')
end
rts = rts/60;

%% runtime by queue
figure(1);clf
hold on

for i = 1:length(viewtypes)
    crts = rts(whvt==i);
    barerr(i,mean(crts),std(crts))
end
set(gca,'XTick',1:length(viewtypes),'XTickLabel',viewtypes)

% for i = 1:length(queues)
%     crts = rts(whq==i);
%     barerr(i,mean(crts),std(crts))
% end
% set(gca,'XTick',1:length(queues),'XTickLabel',queues)

%% runtime over time
% [time,I] = sort(time);
% rts = rts(I);

% figure(2);clf
% plot(time,rts/60,'.')