function panoconvs_patternbits
    totrim = 15;
    outdname = fullfile(mfiledir,'../../figures/patternbits');
    pattdname = fullfile(mfiledir,'../../data/antoinestim/touse');

    if ~exist(outdname,'dir')
        mkdir(outdname);
    end

    d = dir(fullfile(pattdname,'*.png'));
    for i = 1:length(d)
        fprintf('%d/%d:\n',i,length(d))
        
        fn = d(i).name;
        im = rgb2gray(imread(fullfile(pattdname,fn)));
        im = im(totrim+1:end-totrim,:);

        fpre = fullfile(outdname,fn(1:5));
        pattwrite([fpre 'PATT1.png'],im(:,45+(1:90)));
        pattwrite([fpre 'PATT2.png'],im(:,135+(1:90)));
        
        disp('.')
    end
end

function pattwrite(fn,im)
    fprintf('Writing to %s...\n',fn);
    imwrite(im,fn);
end