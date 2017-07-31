clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();

dir_MTSD =[filedir, '/summary/MTSD']; 
sum_dir = [filedir, '/summary'];

%% Number of files to analyse
cd(dir_MTSD);
files = dir('*.csv');
Averages_MTSDall = zeros(1,2);
data_k = zeros(1,2);
Averages_MTSDbinned = zeros(1,4);

for i=1:length(files)
    Name = files(i).name;
    Number = sscanf(Name, '%f');
    Data_name = [num2str(Number),'_summary_MTSD.csv'];
    Data = csvread(Data_name,1,0);
    temp = [Data(:,3), Data(:,5)];
    Data2 = Averages_MTSDall;
    Averages_MTSDall  = [Data2; temp];
end
Averages_MTSDall(1,:) = [];
Averages_MTSDall = sortrows(Averages_MTSDall,1);
cd(sum_dir);
csvwrite('MTSDall.scv', Averages_MTSDall);

for k=0.92:0.005:0.99
    data_k =  Averages_MTSDall( Averages_MTSDall(:,1) >= k-0.005 & Averages_MTSDall(:,1) < k+0.005,:);
    N = int8((k-0.915)/0.005);
    Averages_MTSDbinned(N,1) = k;
    Averages_MTSDbinned(N,2) = 0;
    Averages_MTSDbinned(N,3) = 0;
    if size(data_k,1)>1
        Averages_MTSDbinned(N,2) = mean(data_k(:,2));
        Averages_MTSDbinned(N,3) = sqrt(var(data_k(:,2))/length(data_k(:,2)));
    elseif size(data_k,1)==1
        Averages_MTSDbinned(N,2) = data_k(:,2);
        Averages_MTSDbinned(N,3) = 0;
    end
    Averages_MTSDbinned(N,4) = size(data_k,1);
end

headers = {'Eccentricity', 'MTSD', 'sem','cells'};
csvwrite_with_headers('MTSD_binned.csv',Averages_MTSDbinned,headers);
cd(currdir);
clear variables
clc