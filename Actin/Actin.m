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
act_dir = [filedir, '/actin'];
dens_dir = [filedir, '/cytoskeleton_average'];
MTSD_dir = [filedir, '/cytoskeleton'];
b_dir = [filedir, '/borders'];

%Folder to save information about cells
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];

if exist([filedir,'/summary'],'dir') == 0
    mkdir(filedir,'/summary');
end
sum_dir = [filedir, '/summary'];

if exist([filedir, '/distribution'],'dir') == 0
    mkdir(filedir,'/distribution');
end
dist_dir = [filedir, '/distribution'];



%% Number of files to analyse
cd(dens_dir);
files = dir('*.tif');
cd(currdir);

Averages = zeros(length(files),7);
headers = {'Embryo', 'MT signal area', 'sem','MT bundling','sem',...
    'actin signal area', 'sem','actin signal intensity','sem','Area','sem', ...
    'Eccentricity', 'sem','Direction_cell','sem', 'MT direction','sem', 'direction actin','sem',...
    'Cell number'};

headers2 = {'Cell', 'MT signal area', 'MT  bundling','actin signal area','actin signal intensity',...
    'Area', 'Eccentricity','Direction_cell', ...
    'MT direction', 'direction actin'};

headers3 = {'Embryo','Cell', 'MT signal area', 'MT  bundling','actin signal area','actin signal intensity',...
    'Area', 'Eccentricity','Direction_cell', ...
    'MT direction', 'direction actin'};


Averages_filename = 'Summary.csv';
Averages_filename2 = 'Summary_pulled.csv';
Pulled_data = zeros(1,11);
Dist_filename = 'Distributions.csv';
if exist([filedir,'/summary/cellbycell'],'dir') == 0
    mkdir(filedir,'/summary/cellbycell');
end
result_dir = [filedir, '/summary/cellbycell'];

%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';

all_added_norm = zeros(45,1);
all_added_norm(:,1) = bincenter;


for loop=1:length(files)
    %% reading files
    clear Name Number Actin_file Image_actin Path  Image_borders
    cd(dens_dir);
    Actin_file = [num2str(loop),'.tif'];
    Image = imread(Actin_file);
    cd(MTSD_dir);
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(loop),'/'];
    cd(Path);
    Image_borders = imread('handCorrection.tif');
    [im_x, im_y] = size(Image);
    cd(act_dir);
    Image_actin = imread(Actin_file);
    %% Collect data about cells and boundaries
    borders;
    
    %% Open MTs image, adjust and generate cell masks
    imageprocessing;
    
    %% Cell-by-cell analysis
    
    cellbycelldensity;  
    actindensity;
    
    MTSD;
    vonmises;
    
    actinSD;
    actinvonmises;
    
    %% writing down summarised data
    summaries;
    
    
end
Averages = sortrows(Averages,1);
cd(sum_dir);
csvwrite_with_headers(Averages_filename,Averages,headers);
Pulled_data(1,:) = [];
csvwrite_with_headers(Averages_filename2,Pulled_data,headers3);
csvwrite(Dist_filename, all_added_norm);
cd(currdir);

clc
clear variables
close all
