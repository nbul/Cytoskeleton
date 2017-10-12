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

filedir = uigetdir();
cd(filedir);

average_dir =[filedir, '/cytoskeleton_average'];
b_dir = [filedir, '/borders'];
dens_dir = [filedir, '/summary/kmeans'];
sum_dir = [filedir, '/summary'];
%% Number of files to analyse
cd(average_dir);
files = dir('*.tif');
cd(currdir);

for loop=1:length(files)
    %% reading files
    clear Name Number Actin_file Image_actin Path  Image_borders
    cd(average_dir);
    Name = files(loop).name;
    Number = sscanf(Name, '%f');
    Actin_file = [num2str(Number),'.tif'];
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('tracked_bd.png');
    [im_x, im_y] = size(Image2);
    cd(dens_dir);
    Density_file_name = ['Clustering_',num2str(Number),'.csv'];
    Density_file = csvread(Density_file_name,1,0);
    threshold = Density_file(2);
    %% Collect data about cells and boundaries
    borders;
    DensityCluster;
    
    %% writing down summarised data
    SummariesCluster;
    
end
cd(sum_dir);
Averages_dens = sortrows(Averages_dens,1);
headers_dens = {'Embryo', 'Signal area', 'sem','Density','sem','Bundling','sem',...
    'Uniformity','sem', 'UNAAD','sem','Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'Gaps', 'sem', 'Cell number','Outliers', 'Area', 'sem', 'Eccentricity','sem'};
csvwrite_with_headers('Summary_density_post.csv',Averages_dens,headers_dens);
cd(currdir);

clc
clear variables
close all
