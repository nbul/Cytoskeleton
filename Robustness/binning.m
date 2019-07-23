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

color0 = questdlg(strcat('Which stages'),'Settings','12-15', '15', '15');
if strcmp(color0, '15')
    ecc = 0.92:0.005:0.99;
    int = repmat(0.0025,1,15);
    N = 1:1:15;
else
    ecc = [0.7:0.05:0.95,0.98];
    int = [repmat(0.0025,1,15), 0.0005];
    N = 1:1:7;
end

%% Writing down pulled data MTSD
Averages_MTSDall = zeros(1,2);
Averages_MTSDbinned = zeros(1,5);
cd(MTSD_dir);
for loop=1:length(files)
    %MTSD_file = ['Summary_MTSD',num2str(loop),'.csv'];
    Data = csvread(files(loop).name,1,0);
    temp = [Data(:,5), Data(:,7)];
    Data2 = Averages_MTSDall;
    Averages_MTSDall  = [Data2; temp];
end
Averages_MTSDall(1,:) = [];
Averages_MTSDall = sortrows(Averages_MTSDall,1);
cd(sum_dir);
csvwrite('MTSDall.csv', Averages_MTSDall);

for k=1:1:length(N)
    data_k =  Averages_MTSDall( Averages_MTSDall(:,1) >= ecc(k)-int(k) & Averages_MTSDall(:,1) < ecc(k)+int(k),:);
    Averages_MTSDbinned(N(k),1) = 0;
    Averages_MTSDbinned(N(k),2) = 0;
    Averages_MTSDbinned(N(k),3) = 0;
    Averages_MTSDbinned(N(k),4) = 0;
    if size(data_k,1)>1
        Averages_MTSDbinned(N(k),1) = mean(data_k(:,1));
        Averages_MTSDbinned(N(k),2) = sqrt(var(data_k(:,1))/length(data_k(:,1)));
        Averages_MTSDbinned(N(k),3) = mean(data_k(:,2));
        Averages_MTSDbinned(N(k),4) = sqrt(var(data_k(:,2))/length(data_k(:,2)));
    elseif size(data_k,1)==1
        Averages_MTSDbinned(N(k),1) = data_k(:,1);
        Averages_MTSDbinned(N(k),2) = 0;
        Averages_MTSDbinned(N(k),3) = data_k(:,2);
        Averages_MTSDbinned(N(k),4) = 0;
    end
    Averages_MTSDbinned(N(k),5) = size(data_k,1);
end

headers = {'Eccentricity','sem', 'MTSD', 'sem','cells'};
csvwrite_with_headers('MTSD_binned.csv',Averages_MTSDbinned,headers);

cd(currdir);
clc
clear variables
close all