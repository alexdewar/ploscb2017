jobid = getenv('JOB_ID');

dfile = sprintf('%s/ablate_data%s.mat',mfiledir,jobid);
if ~exist(dfile,'file')
    savedat = true;
    progbaron = true;
    
%     maxinputs = 50; % inc extra generated inputs
    maxcombs = 50;
    nhidden = 10; % hidden units
    ntraincycles = 100;
    ptrain = .4;

    suffix = '_ablate';
    ellfn = sprintf('%s/ellblob%s%s.mat',mfiledir,suffix,jobid);
    if exist(ellfn,'file')
        load(ellfn,'el','az','orients','areas','x','r2_ind','r4_ind','r2ex_ind','r4ex_ind','r24ex_ind','trueimsz');
    else
        genellblobs_ablate;
    end

    t = [orients, areas, el];
    t = bsxfun(@rdivide,bsxfun(@minus,t,min(t)),range(t)); % normalise
%     x = bsxfun(@rdivide,bsxfun(@minus,x,min(x)),range(x)); % normalise

    ioff = prod(trueimsz);
    im_ind = 1:ioff;
    conds = { im_ind; ioff+[r2_ind,r2ex_ind]; ioff+[r4_ind,r4ex_ind]; ...
              ioff+[r2_ind,r4_ind,r24ex_ind] };
    lconds = { 'raw views', 'R2', 'R4d', 'R2+R4d' };
    [ncombs,scores,inputn,inputscores] = deal(cell(length(conds),1));

    for i = 1:length(conds)
        ncombs{i} = NaN(1,length(conds{i}));
        for j = 1:length(conds{i})
            ncombs{i}(j) = nchoosek(length(conds{i}),j);
        end
    end
    
    if progbaron
        startprogbar(10,sum(cellfun(@(p)sum(min(p,maxcombs)),ncombs)));
    end
    for i = 1:length(conds)
        ind = conds{i};

        scores{i} = zeros(1,length(ind));
        [inputscores{i},inputn{i}] = deal(zeros(length(ind)));
        for j = 1:length(ind)
            disp(j)

            if ncombs{i}(j) <= maxcombs
                combs = combnk(1:length(ind),j);
            else
                combs = NaN(maxcombs,j);
                for k = 1:maxcombs
                    combs(k,:) = randperm(length(ind),j);
                end
            end

            tind = randperm(size(t,1));
            traini = tind(1:ptrain*length(tind));
            testi = tind(ptrain*length(tind)+1:end);
            for k = 1:size(combs,1)
                c = combs(k,:);
                net = mlp(j,nhidden,size(t,2),'linear');
                options = zeros(1,18);
                options(14) = ntraincycles;
                net = netopt(net,options,x(traini,ind(c)),t(traini,:),'scg');
                y = mlpfwd(net,x(testi,ind(c)));

                cscore = sum(mean((y-t(testi,:)).^2));
                scores{i}(j) = scores{i}(j) + cscore./size(combs,1);
                inputscores{i}(j,c) = inputscores{i}(j,c) + cscore;
                inputn{i}(j,c) = inputn{i}(j,c)+1;
                if progbaron && progbar
                    return;
                end
            end
        end
        inputscores{i} = inputscores{i}./inputn{i};
    end
    
    if savedat
        save(dfile);
    end
end