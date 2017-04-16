%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and choosing file
data = struct([]);

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
filename = uigetfile('*.csv');
filedata = csvread(filename);
cd(currdir);

filenumber = length(unique(filedata(:,10)));
result = zeros(filenumber, 4);

for i=1:filenumber
    data{i} = filedata(filedata(:,10)==i,:);
    data{i}(:,16) = 1./sqrt(1-data{i}(:,1).*data{i}(:,1));
    data{i}(:,17) = 100*(erf(10./data{i}(:,2)/sqrt(2))-erf(-10./data{i}(:,2)/sqrt(2)))/2;
    result(i,1) = sum(data{i}(:,16))/size(data{i},1);
    result(i,2) = std(data{i}(:,16))/sqrt(size(data{i},1)-1);
    result(i,3) = sum(data{i}(:,17))/size(data{i},1);
    result(i,4) = std(data{i}(:,17))/sqrt(size(data{i},1)-1);
end

summary_filename = ['result_',filename];
headers = {'Aspect ratio','sem', 'Alignment','sem'};
cd(filedir);
csvwrite_with_headers(summary_filename,result,headers);
cd(currdir);
