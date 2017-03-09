r = 123e-3/2;
nim = 500;

load('ofstad_etal_arena.mat');

sq = round(2*sqrt(nim/pi));
[xg,yg] = meshgrid(linspace(-r,r,sq));
selg = hypot(xg,yg)<r;
xg = xg(selg);
yg = yg(selg);

views = NaN(17,90,length(xg));
startprogbar(1,length(xg));
for i = 1:length(xg)
    views(:,:,i) = getalView(xg(i),yg(i),0,0,X,Y,Z);
    if progbar
        return
    end
end

save('ofstad_views.mat','xg','yg','views');