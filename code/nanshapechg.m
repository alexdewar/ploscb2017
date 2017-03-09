nind = [0;find(isnan(X))];

biggest = max(diff(nind));

[X2,Y2,Z2] = deal(NaN(length(nind)-1,biggest));
for i = 1:length(nind)-1
    cn = nind(i)+1:nind(i+1)-1;
    pd = biggest-length(cn);
    X2(i,:) = [X(cn)',X(cn(end))*ones(1,pd)];
    Y2(i,:) = [Y(cn)',Y(cn(end))*ones(1,pd)];
    Z2(i,:) = [Z(cn)',Z(cn(end))*ones(1,pd)];
end