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
%Folders with images
ck_dir =[filedir, '/cytoskeleton'];
dens_dir =[filedir, '/cytoskeleton_average'];
b_dir = [filedir, '/borders'];

%Folder to save information about cells
mkdir(filedir,'/distribution');
dist_dir = [filedir, '/distribution'];

mkdir(filedir,'/images_analysed');
im_dir = [filedir, '/images_analysed'];

mkdir(filedir,'/summary');
sum_dir = [filedir, '/summary'];
%% Number of files to analyse
cd(ck_dir);
files = dir('*.tif');
cd(currdir);


Averages = zeros(length(files),28);
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';

for loop=1:length(files);
    %% reading files
    cd(ck_dir);
    clear Name Number Actin_file Image_actin Path  Image_borders
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
    [im_x, im_y] = size(Image);
    
    %% Collect data about cells and boundaries
    borders;
    
    %% Open MTs image, adjust and generate cell masks
    imageprocessing;
    
    %% Cell-by-cell analysis
    cellbycell;
    cellbycelldensity;
    
    %% Von-Mises fit
    vonmises_fit_dist_sum;
    
    %% writing down summarised data
    summaries;
    
end
Averages = sortrows(Averages,1);
summary_filename = 'Summary_all.csv';
headers = {'Cell', 'Density', 'sem', 'SD', 'sem', 'Direction_cytoskeleton','sem', 'Area','sem', ...
    'Eccentricity','sem', 'Direction_cell','sem', 'DEV','sem', 'Signal Area','sem',...
    'Aspect ratio','sem', 'Alignment','sem', 'Sparseness','sem','Bundling','sem',...
    'Uniformity','sem','Cell number'};
cd(sum_dir);
csvwrite_with_headers(summary_filename,Averages,headers);
cd(currdir);

clc
clear variables
close all
