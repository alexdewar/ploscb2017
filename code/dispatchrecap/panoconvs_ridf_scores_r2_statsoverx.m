clear

load('panoconv_nothresh.mat')
x = eh_dcp;
y = diffr2(~isnan(x));
x = x(~isnan(x));
[dat,xp] = StatsOverX(x,y,linspace(min(x),max(x),10))