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


%% Determining mode and method
method = 1;
path = 1;

usedefault = questdlg(strcat('Which type of analysis?'),'Settings','MTSD','Density','MTSD');
if strcmp(usedefault, 'Density')
    path = 0;
end
if path == 0
    usedefault = questdlg(strcat('Which threshold methos?'),'Settings','kmeans','Otsu','2DOtsu','kmeans');
    if strcmp(usedefault, 'Otsu')
        method = 0;
    elseif strcmp(usedefault, '2DOtsu')
        method = 2;
    elseif strcmp(usedefault, 'kmeans')
        method = 3;
    end
end

%% Number of files to analyse
cd(ck_dir);
files = dir('*.tif');
cd(currdir);

if path ==1
    Averages = zeros(length(files),17);
    headers = {'Embryo', 'Area','sem',  'Eccentricity','sem', 'Direction_cell','sem',...
        'SD', 'sem', 'Direction_cytoskeleton','sem', 'Aspect ratio','sem', 'Alignment','sem','Cell number','Outliers'};
    Averages_filename = 'Summary_MTSD.csv';
    if exist([filedir,'/summary/MTSD'],'dir') == 0
        mkdir(filedir,'/summary/MTSD');
    end
    result_dir = [filedir, '/summary/MTSD'];
else
    Averages = zeros(length(files),17);
    headers = {'Embryo', 'Signal area', 'sem','Density','sem','Bundling','sem', 'Uniformity','sem', ...
        'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem','Cell number','Outliers'};
    if method == 1
        Averages_filename = 'Summary_edges.csv';
         if exist([filedir,'/summary/edges'],'dir') == 0
             mkdir(filedir,'/summary/edges');
         end
         result_dir = [filedir, '/summary/edges'];
    elseif method == 0
         Averages_filename = 'Summary_otsu.csv';
         if exist([filedir,'/summary/otsu'],'dir') == 0
             mkdir(filedir,'/summary/otsu');
         end
         result_dir = [filedir, '/summary/otsu'];
    elseif method == 2
        Averages_filename = 'Summary_2Dotsu.csv';
        if exist([filedir,'/summary/2Dotsu'],'dir') == 0
            mkdir(filedir,'/summary/2Dotsu');
        end
        result_dir = [filedir, '/summary/2Dotsu'];
    elseif method == 3
        Averages_filename = 'Summary_kmeans.csv';
        if exist([filedir,'/summary/kmeans'],'dir') == 0
            mkdir(filedir,'/summary/kmeans');
        end
        result_dir = [filedir, '/summary/kmeans'];
    end
end
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
    if path == 1
        cellbycell;
        %% Von-Mises fit
        vonmises_fit_dist_sum;
    else
        cellbycelldensity;
    end
    
    
    %% writing down summarised data
    summaries;
    
end
Averages = sortrows(Averages,1);
% summary_filename = 'Summary_all.csv';
% headers = {'Cell', 'Density', 'sem', 'SD', 'sem', 'Direction_cytoskeleton','sem', 'Area','sem', ...
%     'Eccentricity','sem', 'Direction_cell','sem', 'DEV','sem', 'Signal Area','sem',...
%     'Aspect ratio','sem', 'Alignment','sem', 'Sparseness','sem','Bundling','sem',...
%     'Uniformity','sem','Cell number'};
cd(sum_dir);
csvwrite_with_headers(Averages_filename,Averages,headers);
cd(currdir);

clc
clear variables
close all
