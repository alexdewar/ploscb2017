function [nx,ny] = appendshapes(varargin)
midi = floor(length(varargin)/2);
x = varargin(1:midi);
y = varargin(midi+1:end);

len = max(cellfun(@(a)size(a,1),x));
wid = cellfun(@(a)size(a,2),x);
cwid = cumsum(wid);
[nx,ny] = deal(NaN(len,sum(wid)));
for i = 1:numel(x)
    cur = cwid(i):cwid(i)+wid(i)-1;
    ind = 1:length(x{i});
    nx(ind,cur) = x{i};
    ny(ind,cur) = y{i};
    
    ind2 = length(x{i})+1:len;
    nx(ind2,cur) = x{i}(end,:);
    ny(ind2,cur) = y{i}(end,:);
end