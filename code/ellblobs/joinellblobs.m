nfile = 10;

d = dir([mfiledir '/ellblob*_*.mat']);

load([mfiledir d(1).name]);
x_im