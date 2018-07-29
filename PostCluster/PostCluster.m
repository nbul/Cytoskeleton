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

folders;

%% Number of files to analyse
cd(dens_dir);
files = dir('*.tif');
cd(currdir);
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';
stretch_factor = 0.5;
for loop=1:length(files)
    %% reading files
    cd(MTSD_dir);
    Name = files(loop).name;
    Number = sscanf(Name, '%f');
    Actin_file = [num2str(Number),'.tif'];
    Image = imread(Actin_file);
    cd(dens_dir);
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('tracked_bd.png');
    [im_x, im_y] = size(Image2);
    cd(th_dir);
    Density_file_name = ['Clustering_',num2str(Number),'.csv'];
    Density_file = csvread(Density_file_name,1,0);
    threshold = Density_file(2);
    %% Collect data about cells and boundaries
    borders;
    images;
    
    MTSD;
    vonmises;    
    Density;
    
    %% writing down summarised data
    Summary;
    
end
cd(sum_dir);
Averages_dens = sortrows(Averages_dens,1);
Averages_MTSD = sortrows(Averages_MTSD,1);
csvwrite_with_headers('Summary_density.csv',Averages_dens,headers_dens_all);
csvwrite_with_headers('Summary_MTSD.csv',Averages_MTSD,headers_MTSD_all);
cd(currdir);

clc
clear variables
close all
