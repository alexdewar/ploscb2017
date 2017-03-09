clear

njob=10;
nvar=1000;
nwave=2;

n = njob*nvar;

[thoffs,scales,mjas,mnas,rat] = deal(NaN(n,1));
d = dir([mfiledir '/ellfp*.mat']);

for i = 1:length(d)
    load([mfiledir '/' d(i).name]);
    
    ind = starti+(0:nvar-1);
    thoffs(ind) = thoff;
    scales(ind) = scale;
    rat(ind) = amp(:,1)./amp(:,2);
    mjas(ind) = mja;
    mnas(ind) = mna;
end

cpar = rat;
figure(1);clf
[cpar,I] = sort(cpar);
plot(cpar,mjas(I))