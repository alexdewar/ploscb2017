function stats_R2vsSPP

showprogbar = true;
usenothreshkerns = true;

showrng = [45 135 225];
im_size = [120 360];

dname = '../data/patternstim';

if usenothreshkerns
    threshstr = '_nothresh';
else
    threshstr = '';
end

load(['vf_kernels' threshstr],'vf_avkernels_r2');
kerns = kstr2k(vf_avkernels_r2);
if usenothreshkerns
    for i = 1:size(kerns,3)
        ck = kerns(:,:,i);
        
        pos = ck > 0;
        ck(pos) = ck(pos)./sum(ck(pos));
        neg = ck < 0;
        ck(neg) = -ck(neg)./sum(ck(neg));
        
        kerns(:,:,i) = ck;
    end
end

fulldname = fullfile(mfiledir,dname,'touse');
d = [dir(fullfile(fulldname,'*.jpg'));dir(fullfile(fulldname,'*.png'))];
fname = sort({d.name});

crng = round(showrng*im_size(2)/360);
[straightdiffr2,stdr2,eh_scp] = deal(NaN(length(fname),1));

datafn = fullfile(mfiledir, 'panoconv_R2_straightdiff.mat');
if exist(datafn,'file')
    load(datafn);
else
    if showprogbar
        startprogbar(1,length(fname));
    end
    
    ehf = fopen([mfiledir,'/dcp.txt'],'r');
    ehcell = textscan(ehf,'%d:%s\n');
    fclose(ehf);
    eh_fig = ehcell{1};
    
    for i = 1:length(fname)
        im = imread(fullfile(fulldname,fname{i}));
        if size (im,3)>1;
            im = rgb2gray(im);
        end
        im = im2double(im);
        
        eh_scp(i) = str2double(fname{i}(13))*0.1*str2double(fname{i}(15:18))/str2double(fname{i}(9:11));
        
        if ~isnan(eh_scp(i))
            [straightdiffr2(i),stdr2(i)] = meandiff(im,kerns);
        end
        
        if showprogbar && progbar
            return;
        end
    end
    save(datafn,'straightdiffr2','stdr2','eh_fig','eh_scp');
end

%% do stats
valid = ~isnan(eh_scp);
X = straightdiffr2(valid);
Y = eh_scp(valid);

fprintf('N(valid) = %d\nN(total) = %d\n', sum(valid), length(eh_scp));

[rho,pval] = corr(X,Y,'type','Spearman')

dump2base(true)

function [dr2,stdr2]=meandiff(im,kerns)
ksz = size(kerns);
rim = imresize(im,[ksz(1),ceil(360*ksz(2)/270)]);
xoff = (size(rim,2)-ksz(2))/2;
acts = NaN(ksz(3),size(rim,2));
for i = 1:size(rim,2)
    crim = circshift(rim,[0 i-1-ceil(size(rim,2)/2)]);
    
    acts(:,i) = shiftdim(sum(sum(bsxfun(@times,crim(:,xoff+(1:ksz(2))),kerns)),2));
end

acts0 = acts(:,1+size(rim,2)/2);
rr2 = mean(bsxfun(@minus,acts,acts0));
dr2 = mean(rr2);

stdr2 = std(rr2)/sqrt(length(rr2));