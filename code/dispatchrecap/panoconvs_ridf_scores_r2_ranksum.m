clear

load('panoconv_nothresh.mat','diffr2','eh_dcpsig')
testable = ~isnan(eh_dcpsig);
sigs = testable & eh_dcpsig > 0;
fprintf('N(sig) = %d\n', sum(sigs));
notsigs = testable & eh_dcpsig == 0;
fprintf('N(notsig) = %d\n', sum(notsigs));

% h = kstest(diffr2(testable))
% [hsw, p, swstat] = swtest(diffr2(testable))
% [p,h,stats] = ranksum(diffr2(sigs),diffr2(notsigs))
