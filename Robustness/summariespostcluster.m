%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
sum_dir = [filedir, '/summary'];
MTSD_dir = [filedir, '/summary/cellbycell'];

%% Writing down pulled data MTSD
Averages_MTSDall = zeros(1,2);
data_k = zeros(1,2);
Averages_MTSDbinned = zeros(1,5);
cd(MTSD_dir);
files = dir('*.csv');
for loop=1:length(files)
    Data = csvread(files(loop).name,1,0);
    temp = [Data(:,4), Data(:,6)];
    Data2 = Averages_MTSDall;
    Averages_MTSDall  = [Data2; temp];
end
Averages_MTSDall(1,:) = [];
Averages_MTSDall = sortrows(Averages_MTSDall,1);
cd(sum_dir);
csvwrite('MTSDall.csv', Averages_MTSDall);

for k=0.92:0.005:0.99
    data_k =  Averages_MTSDall( Averages_MTSDall(:,1) >= k-0.0025 & Averages_MTSDall(:,1) < k+0.0025,:);
    N = int8((k-0.915)/0.005);
    Averages_MTSDbinned(N,1) = 0;
    Averages_MTSDbinned(N,2) = 0;
    Averages_MTSDbinned(N,3) = 0;
    Averages_MTSDbinned(N,4) = 0;
    if size(data_k,1)>1
        Averages_MTSDbinned(N,1) = mean(data_k(:,1));
        Averages_MTSDbinned(N,2) = sqrt(var(data_k(:,1))/length(data_k(:,1)));
        Averages_MTSDbinned(N,3) = mean(data_k(:,2));
        Averages_MTSDbinned(N,4) = sqrt(var(data_k(:,2))/length(data_k(:,2)));
    elseif size(data_k,1)==1
        Averages_MTSDbinned(N,1) = data_k(:,1);
        Averages_MTSDbinned(N,2) = 0;
        Averages_MTSDbinned(N,3) = data_k(:,2);
        Averages_MTSDbinned(N,4) = 0;
    end
    Averages_MTSDbinned(N,5) = size(data_k,1);
end

headers = {'Eccentricity','sem', 'MTSD', 'sem','cells'};
csvwrite_with_headers('MTSD_binned.csv',Averages_MTSDbinned,headers);

cd(currdir);
clc
clear variables
close all