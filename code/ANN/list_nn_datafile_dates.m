function list_nn_datafile_dates

dname = 'drosodata/ANNs';
d = dir(fullfile(dname,'figpreprocess_*.mat'));

for i = 1:length(d)
    load(fullfile(dname,d(i).name),'metadata','fname_ellblob')
    fprintf('\n%s:\n',d(i).name)
    disp(fname_ellblob)
    disp(metadata)
end