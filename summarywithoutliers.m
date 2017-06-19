clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
ck_dir =[filedir, '/summary/outliers']; 

%% Number of files to analyse
cd(ck_dir);
files = dir('*.csv');
cd(currdir);
Averages = zeros(length(files),28);
for i=1:numel(files)
    cd(ck_dir);
    Name = files(i).name;
    Number = sscanf(Name, '%f');
    Data_name = [num2str(Number),'_summary.csv'];
    Data = csvread(Data_name,1,0);
    
    Averages(i,1) = Number;
    % Density
    Averages(i,2) = nanmean(Data(:,2));
    Averages(i,3) = sqrt(var(Data(:,2))/length(Data(:,2)));
    %SD
    Averages(i,4) = mean(Data(:,3));
    Averages(i,5) = sqrt(var(Data(:,3))/length(Data(:,3)));
    %Direction cytoskeleton
    Averages(i,6) = mean(Data(:,4));
    Averages(i,7) = sqrt(var(Data(:,4))/length(Data(:,4)));
    %Area
    Averages(i,8) = mean(Data(:,5));
    Averages(i,9) = sqrt(var(Data(:,5))/length(Data(:,5)));
    % Eccentricity
    Averages(i,10) = mean(Data(:,6));
    Averages(i,11) = sqrt(var(Data(:,6))/length(Data(:,6)));
    % Direction_cell
    Averages(i,12) = mean(Data(:,7));
    Averages(i,13) = sqrt(var(Data(:,7))/length(Data(:,7)));
    %DEV
    Averages(i,14) = mean(Data(:,8));
    Averages(i,15) = sqrt(var(Data(:,8))/length(Data(:,8)));
    %Signal area
    Averages(i,16) = mean(Data(:,9));
    Averages(i,17) = sqrt(var(Data(:,9))/length(Data(:,9)));
    %Aspect ratio
    Averages(i,18) = mean(Data(:,10));
    Averages(i,19) = sqrt(var(Data(:,10))/length(Data(:,10)));
    %Alignment
    Averages(i,20) = mean(Data(:,11));
    Averages(i,21) = sqrt(var(Data(:,11))/length(Data(:,11)));
    %Sparseness
    Averages(i,22) = nanmean(Data(:,12));
    Averages(i,23) = sqrt(var(Data(~isnan(Data(:,12)), 12))/...
        length(Data(~isnan(Data(:,12)), 12)));
    %Bundling
    Averages(i,24) = mean(Data(:,13));
    Averages(i,25) = sqrt(var(Data(:,13))/length(Data(:,13)));
    %Bundling
    Averages(i,26) = mean(Data(:,14));
    Averages(i,27) = sqrt(var(Data(:,14))/length(Data(:,13)));
    
    %number of cells
    Averages(i,28) = length(Data(:,11));
end
Averages = sortrows(Averages,1);
summary_filename = 'Summary_all_clean.csv';
headers = {'Embryo', 'Density', 'sem', 'SD', 'sem', 'Direction_cytoskeleton','sem', 'Area','sem', ...
    'Eccentricity','sem', 'Direction_cell','sem', 'DEV','sem', 'Signal Area','sem',...
    'Aspect ratio','sem', 'Alignment','sem', 'Sparseness','sem','Bundling' ,'sem',...
    'Uniformity','sem','Cell number'};
csvwrite_with_headers(summary_filename,Averages,headers);
cd(currdir);

clc
clear variables
close all