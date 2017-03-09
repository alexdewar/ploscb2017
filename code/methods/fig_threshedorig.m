whichkernsforfig;

%% show stuff
figure(1);clf
alsubplot(3,2,1,1)
imshow(kernfn);

alsubplot(1,2)
showkernel(kern(1));
axis equal off

alsubplot(2,1)
hold on
ac = mean([kern(1).cent;kern(2).cent]);
showkernels(kern,[],kalpha);
plot(ac(1),ac(2),'y+')
hold off
axis equal off

alsubplot(2,2)
ckern = centerkernson(kern,ac);
showkernels(ckern,[],kalpha);
axis equal off

alsubplot(3,1)
acs = mean(cell2mat({kerns.cent}'));
ckerns = centerkernson(kerns,acs);
showkernels(ckerns,[],kalpha/2);
axis equal off

alsubplot(3,2)
showkernel(centerkernson(avkern,ac));
axis equal off

%savefig('averaging',[8 8],'png');