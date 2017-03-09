dname = [mfiledir '/../../data/antoinestim/touse'];
outdname = '/home/alex/sync/neuroscience day 2015 talk/patterns';

d = dir([dname '/*.png']);

for i = 1:length(d)
    im = imread(fullfile(dname,d(i).name));
    imwrite(im(:,45+(1:180),:),fullfile(outdname,d(i).name));
end