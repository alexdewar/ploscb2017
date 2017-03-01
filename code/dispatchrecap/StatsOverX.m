% function[dat,xp]=StatsOverX(x,y,xdiv)
%
% this function generates various stats based on y after binning the data
% based on x 
%
% xdiv allows you to specify the edges of the bins (if a vector), or a
% number of bins (if a number) and defaults to 20
%
% it returns a structure dat which has various fields and xp which is the
% centre of the bins used (good for plotting
%
% usage: to plot mean and sd:
% errorbar(xp,[dat.me],[dat.sd])
%
% to plot median and 25th and 75th
% errorbar(xp,[dat.med],[dat.med]-[dat.p25],[dat.p75]-[dat.med])
% 
% to plot circular mean and sd:
%  errorbar(xp,[dat.meang],[dat.angsd])
 
function[dat,xp]=StatsOverX(x,y,xdiv)

if(nargin<3) xdiv=20; end;
if(length(xdiv)==1)
    minx=min(x);
    maxx=max(x);
    wid=(maxx-minx)/xdiv;
    xdiv=minx:wid:maxx;
end

for i=1:(length(xdiv)-1)
    is=(x>=xdiv(i))&(x<xdiv(i+1));
    xp(i)=0.5*(xdiv(i)+xdiv(i+1));
    d=y(is);
    dat(i).med=median(d);    
    dat(i).me=mean(d);
    [mea,vl]=circ_mean(d);
    dat(i).meang=mea;
    dat(i).meangL=vl;
    if(length(d>1))
        dat(i).medang=circ_median(d');
    else
        dat(i).medang=circ_median(d);
    end
    dat(i).angsd=sqrt(2*(1-vl));
    dat(i).sd=std(d);
    dat(i).p25=prctile(d,25);
    dat(i).p75=prctile(d,75);
    dat(i).n=length(d);
    dat(i).x=d;
    dat(i).xdiv=xdiv;
end    