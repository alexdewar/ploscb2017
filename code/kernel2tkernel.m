function tkerns=kernel2tkernel(kerns,cents)
    if nargin<2
        cents = cell2mat({kerns.cent}');
        if isfield(kerns,'isleft')
            lefts = cell2mat({kerns.isleft});
        end
        
        kerns = {kerns.k};
    else
        ckern = kerns;
    end
    
    tkerns(size(cents,1)) = struct('k',[],'cent',[],'isleft',[],'oldcent',[],'im_size',[]);
    for i = 1:length(tkerns)
        if iscell(kerns)
            ckern = kerns{i};
        end
        tkerns(i).oldcent = cents(i,:);
        [tkerns(i).k,xtrim,ytrim] = bwtrim(ckern);
        tkerns(i).cent = cents(i,:)-[xtrim(1) ytrim(1)]+1;
        tkerns(i).im_size = size(ckern);
        if exist('lefts','var')
            tkerns(i).isleft = lefts(i);
        end
    end
end