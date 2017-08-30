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
dens_dir =[currdir, '/cytoskeleton_average'];
b_dir = [currdir, '/borders'];

mkdir(currdir,'/summary');
sum_dir = [currdir, '/summary'];
method = 3;

%% Number of files to analyse
cd(dens_dir);
files = dir('*.tif');
cd(currdir);
Averages_filename = 'Summary_kmeans.csv';
mkdir(currdir,'/summary/kmeans');
result_dir = [currdir, '/summary/kmeans'];
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';
warning('off','stats:kmeans:FailedToConvergeRep');
for loop=1:length(files)
    %% reading files
    clear Name Number Actin_file Image_actin Path  Image_borders
    cd(dens_dir);
    Name = files(loop).name;
    Number = sscanf(Name, '%f');
    Actin_file = [num2str(Number),'.tif'];
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('tracked_bd.png');
    [im_x, im_y] = size(Image);
    
    %% Collect data about cells and boundaries
    borders;
    cellbycelldensity;
    
    %% writing down summarised data
    summaries_cluster;
    
end
Averages = sortrows(Averages,1);
% summary_filename = 'Summary_all.csv';
headers = {'Embryo', 'Signal area', 'sem','Density','sem','Bundling','sem', 'Uniformity','sem', ...
        'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem','Cell number'};
cd(sum_dir);
csvwrite_with_headers(Averages_filename,Averages,headers);
cd(currdir);

clc
clear variables
close all
