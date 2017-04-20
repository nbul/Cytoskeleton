%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and choosing file


currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
filename = uigetfile('*.csv');
filedata = csvread(filename);
cd(currdir);
data = 100*(erf(10./filedata/sqrt(2))-erf(-10./filedata/sqrt(2)))/2;


result(:,1) = [1.400280084; 1.511857892; 1.666666667; 1.898315992; 2.294157339;3.202563076;5.025189076];
result(:,2) = sum(data,1)'/size(data,1)';
result(:,3) = std(data,1)'/sqrt(size(data,1)-1)';

summary_filename = ['result_',filename];
headers = {'Aspect ratio','Alignment','sem'};
cd(filedir);
csvwrite_with_headers(summary_filename,result,headers);
cd(currdir);
