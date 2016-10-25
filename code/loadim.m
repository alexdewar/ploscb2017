function im=loadim(imfn)
[pth,fn,ext] = fileparts(imfn);
if strcmpi(ext,'.mat')
    vars = who('-file',imfn);
    if numel(vars)~=1
        error('needs to be only one var in file');
    else
        load(imfn,vars{1});
        eval(['im=' vars{1} ';']);
    end
else
    im = imread(imfn);
end