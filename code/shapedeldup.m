function [nx,ny] = shapedeldup(x,y)
isrep = (circdiff(x)==0 & circdiff(y)==0);
nrep = sum(isrep,1);
minrep = min(nrep);
nlst = nrep-minrep;
len = size(x,1)-minrep;
[nx,ny] = deal(NaN(len,size(x,2)));
for i = 1:size(x,2)
    if nlst(i) > 0
        lastrep = find(isrep(:,i),nlst(i),'last');
        isrep(lastrep,i) = false;
    end
    nx(:,i) = x(~isrep(:,i),i);
    ny(:,i) = y(~isrep(:,i),i);
end