load worldtrees2

ntest = 10;

xs = (range(X(:))*rand(ntest,1)-min(X(:)))/2;
ys = (range(Y(:))*rand(ntest,1)-min(Y(:))/2);
zs = (max(Z(:))*rand(ntest,1))/2;
ths = 2*pi*rand(ntest,1);

for i = 1:ntest
    v = getView(xs(i),ys(i),zs(i),ths(i),X,Y,Z,true,false,true);
    
    figure(i);clf
    hold on
    subplot(2,1,1)
    imshow(getviewfast(xs(i),ys(i),zs(i),ths(i),X,Y,Z,true));

    subplot(2,1,2);
    imshow(v);
end