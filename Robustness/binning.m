clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

dist_dir = [filedir, '/distribution'];
sum_dir = [filedir, '/summary'];
MTSD_dir = [filedir, '/summary/cellbycell'];

%% Writing down summarised data MTSD and density embryo by embryo
cd(MTSD_dir);
files = dir('*.csv');


%% Writing down pulled data MTSD
Averages_MTSDall = zeros(1,2);
data_k = zeros(1,2);
Averages_MTSDbinned = zeros(1,5);
cd(MTSD_dir);
for loop=1:length(files)
    %MTSD_file = ['Summary_MTSD',num2str(loop),'.csv'];
    Data = csvread(files(loop).name,1,0);
    temp = [Data(:,4), Data(:,6)];
    Data2 = Averages_MTSDall;
    Averages_MTSDall  = [Data2; temp];
end
Averages_MTSDall(1,:) = [];
Averages_MTSDall = sortrows(Averages_MTSDall,1);
cd(sum_dir);
csvwrite('MTSDall.csv', Averages_MTSDall);

for k=0.7:0.05:0.95
    data_k =  Averages_MTSDall( Averages_MTSDall(:,1) >= k-0.025 & Averages_MTSDall(:,1) < k+0.025,:);
    N = int8((k-0.65)/0.05);
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

data_k =  Averages_MTSDall( Averages_MTSDall(:,1) >= 0.97,:);

Averages_MTSDbinned(7,1) = mean(data_k(:,1));
Averages_MTSDbinned(7,2) = sqrt(var(data_k(:,1))/length(data_k(:,1)));
Averages_MTSDbinned(7,3) = mean(data_k(:,2));
Averages_MTSDbinned(7,4) = sqrt(var(data_k(:,2))/length(data_k(:,2)));
Averages_MTSDbinned(7,5) = size(data_k,1);
headers = {'Eccentricity','sem', 'MTSD', 'sem','cells'};
csvwrite_with_headers('MTSD_binned.csv',Averages_MTSDbinned,headers);

cd(currdir);
clc
clear variables
close all