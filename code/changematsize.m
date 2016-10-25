function x=changematsize(x,sz)
chg = sz-size(x);

if chg(1) >= 0
    x = [x;zeros(chg(1),size(x,2))];
else
    x = x(1:sz(1),:);
end

if chg(2) >= 0
    x = [x,zeros(sz(1),chg(2))];
else
    x = x(:,1:sz(2));
end