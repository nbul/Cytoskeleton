%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

%% Parameters
bin_size = 10;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;

cd(filedir);
datafile = csvread('data_angles.csv');
[N, bins] = histc(datafile(:,1),bincenter);

m_added_norm = [bincenter', N];

 vonmises_fit_dist_sum;