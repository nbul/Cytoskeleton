% Script for Automatic Analysis of N embryos
% The number of filesets should be entered for each experiment
% Adapted to use Sobel 5x5 to find direction
% Check analysis of Microtubule density

clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
average_dir =[currdir, '/cytoskeleton_average'];

mkdir(currdir,'/summary');
sum_dir = [currdir, '/summary'];


%% Number of files to analyse
cd(average_dir);
files = dir('*.tif');
cd(currdir);
mkdir(currdir,'/summary/threshold');
dens_dir = [currdir, '/summary/threshold'];

%% Parameters
warning('off','stats:kmeans:FailedToConvergeRep');
for loop=1:length(files)
    %% reading files
    clear Name Number Actin_file Image_actin Path  Image_borders
    cd(average_dir);
    Name = files(loop).name;
    Number = sscanf(Name, '%f');
    Actin_file = [num2str(Number),'.tif'];
    Image2 = imread(Actin_file);
    [im_x, im_y] = size(Image2);
    
    %% Threshold
    DensityCluster;  

    
end
cd(currdir);

clc
clear variables
close all
