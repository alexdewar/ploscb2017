function savefig(ofname,sz,ext,type)
% function savefig(fname,sz,ext,type)
global verbose
if nargin == 3
    type = ext;
end
% issvg = false;
if nargin < 3
    type = 'pdf';
    ext = 'pdf';
% elseif strcmp(type,'svg')
%     issvg = true;
%     type = 'pdf';
%     ext = 'pdf';
end
ext = ['.' ext];
if nargin < 2
    sz = [16 7];
end
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',sz);
set(gcf,'PaperPosition',[0 0 sz]);
% set(gca,'FontSize',8)
%set(gca, 'Position', get(gca, 'OuterPosition') - ...
%   get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
%set(gca,'Position',[0 0 1 1]);

% if length(fname) >= length(ext) && strcmpi(fname((end-length(ext)+1):end),ext)
%     fname = fname((end-length(ext)+1):end);
% % end
% ofname = [figdir filesep ofname];
cnt = 1;
while true
    fname = sprintf('%s (%04d)%s',ofname,cnt,ext);
    if ~exist(fname,'file')
        break;
    end
    cnt = cnt+1;
end
fprintf('Saving figure to %s...\n',fname);
print(['-d' type],'-r600',fname);

% if issvg % do conversion
%     sfname = [figdir filesep ofname];
%     cnt = 2;
%     while exist([sfname '.svg'],'file')
%         sfname = sprintf('%s (%d)',ofname,cnt);
%         cnt = cnt+1;
%     end
%     sfname = [sfname '.svg'];
%     fprintf('Converted file name: %s\n',sfname);
%     system(sprintf('inkscape -l "%s" "%s%s"',sfname,fname,ext));
% end