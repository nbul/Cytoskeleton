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
if exist([filedir, '/distribution'],'dir') == 0
    mkdir(filedir,'/distribution');
end
dist_dir = [filedir, '/distribution'];

if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];

if exist([filedir,'/summary'],'dir') == 0
    mkdir(filedir,'/summary');
end
sum_dir = [filedir, '/summary'];


%% Number of files to analyse
cd(ck_dir);
files = dir('*.tif');
cd(currdir);


Averages_MTSD = zeros(length(files),17);
headers_MTSD = {'Embryo', 'Area','sem',  'Eccentricity','sem', 'Direction_cell','sem',...
    'SD', 'sem', 'Direction_cytoskeleton','sem', 'Aspect ratio','sem', 'Alignment','sem','Cell number','Outliers'};
Averages_filename_MTSD = 'Summary_MTSD.csv';
if exist([filedir,'/summary/MTSD'],'dir') == 0
    mkdir(filedir,'/summary/MTSD');
end
result_dir_MTSD = [filedir, '/summary/MTSD'];

Averages_Dens = zeros(length(files),15);
headers_Dens = {'Embryo', 'Uniformity','sem', 'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'Sdr','sem', 'Sdq','sem','Cell number','Outliers'};
Averages_filename_Dens = 'Summary_Dens.csv';
if exist([filedir,'/summary/Dens'],'dir') == 0
    mkdir(filedir,'/summary/Dens');
end
result_dir_Dens = [filedir, '/summary/Dens'];

%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';

for loop=1:length(files)
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
    %% Von-Mises fit
    vonmises_fit_dist_sum;
    
    cellbycelldensity;
    
    
    %% writing down summarised data
    summaries;
    
end
Averages_MTSD = sortrows(Averages_MTSD,1);
Averages_Dens = sortrows(Averages_Dens,1);
cd(sum_dir);
csvwrite_with_headers(Averages_filename_MTSD,Averages_MTSD,headers_MTSD);
csvwrite_with_headers(Averages_filename_Dens,Averages_Dens,headers_Dens);
cd(currdir);

clc
clear variables
close all
