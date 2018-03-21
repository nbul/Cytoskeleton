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
max_dir =[currdir, '/cytoskeleton'];
average_dir =[currdir, '/cytoskeleton_average'];
b_dir = [currdir, '/borders'];

mkdir(currdir,'/summary');
sum_dir = [currdir, '/summary'];
mkdir(currdir,'/distribution');
dist_dir = [currdir, '/distribution'];

%% Number of files to analyse
cd(max_dir);
files = dir('*.tif');
cd(currdir);
mkdir(currdir,'/summary/kmeans');
dens_dir = [currdir, '/summary/kmeans'];
mkdir(currdir,'/summary/MTSD');
MTSD_dir = [currdir, '/summary/MTSD'];
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
    cd(max_dir);
    Name = files(loop).name;
    Number = sscanf(Name, '%f');
    Actin_file = [num2str(Number),'.tif'];
    Image = imread(Actin_file);
    cd(average_dir);
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('handCorrection.tif');
    [im_x, im_y] = size(Image2);
    
    %% Collect data about cells and boundaries
    borders;
    MTSDcluster;
    vonmises_fit_dist_sum;
    DensityCluster;  
    
    %% writing down summarised data
    SummariesCluster;
    
end
cd(currdir);

clc
clear variables
close all
