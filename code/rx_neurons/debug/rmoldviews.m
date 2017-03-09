curdir = pwd;

cd ~/curdrosproj/data/rx_neurons/views

d = dir('*.mat');

for i = 1:length(d)
    if varsinmatfile(d(i).name,'kviews')
        load(d(i).name,'kviews');
        mns = min(kviews);
        mxs = max(kviews);
        if ~all(mns(:)==0) || ~all(mxs(:)==1)
            disp(d(i).name)
        end
    else
        warning('no kviews found in %s',d(i).name);
    end
end

cd(curdir)