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


Averages = zeros(length(files),23);
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;

for loop=1:length(files);
    %% reading files
    cd(ck_dir);
    clear Name Number Actin_file Image_actin Path  Image_borders
    Name = files(loop).name;
    Number = sscanf(Name, '%f');
    Actin_file = [num2str(Number),'.tif'];
    Image_actin = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('tracked_bd.png');
    
    %% Running analysis script
    actin_analysis_v1;
    
    %% writing down distributions in each image and analysed image
    cd(dist_dir);
    gradient_filename = [num2str(Number),'_distribution.csv'];
    csvwrite(gradient_filename, m_added_norm);
    
    cd(im_dir);
    image_filename = [num2str(Number),'_analysed_image.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    
    %% Von-Mises fit
    vonmises_fit_dist_sum;
    
    %% Writing down summarised data
    
    summary = zeros(length(SD),12);
    for counter2 = 1:length(SD)
        summary(counter2,1) = counter2;
    end
    summary(:,2) = mts_density';
    summary(:,3) = SD';
    summary(:,4) = mu';
    summary(:,5) = cell_data(:,2);
    summary(:,6) = cell_data(:,3);
    summary(:,7) = cell_data(:,4);
    for counter2 = 1:length(SD)
        if (abs(summary(counter2, 4)-summary(counter2, 7)) >= 90)
            summary(counter2, 8) = 180 - abs(summary(counter2, 4)-summary(counter2, 7));
        else
            summary(counter2, 8) = abs(summary(counter2, 4)-summary(counter2, 7));
        end
    end
    summary(:,9) = mts_area';
    summary(:, 10) = 1./sqrt(1-cell_data(:,3).*cell_data(:,3));
    summary(:,11) = 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2;
    summary(:,12) = Spars';
    
    summary_filename = [num2str(Number),'_summary.csv'];
    headers = {'Cell', 'Density', 'SD', 'Direction_cytoskeleton','Area', 'Eccentricity', 'Dorection_cell','DEV', 'Signal Area','Aspect ratio','Alignment', 'Sparseness'};
    cd(sum_dir);
    csvwrite_with_headers(summary_filename,summary,headers);
    
    
    Averages(loop,1) = loop;
    % Density
    Averages(loop,2) = nanmean(mts_density);
    Averages(loop,3) = var(mts_density(~isnan(mts_density)))/sqrt(length(mts_density(~isnan(mts_density)))-1);
    %SD
    Averages(loop,4) = mean(SD);
    Averages(loop,5) = var(SD)/sqrt(length(SD)-1);
    %Direction cytoskeleton
    Averages(loop,6) = mean(mu);
    Averages(loop,7) = var(mu)/sqrt(length(SD)-1);
    %Area
    Averages(loop,8) = mean(cell_data(:,2));
    Averages(loop,9) = var(cell_data(:,2))/sqrt(length(SD)-1);
    % Eccentricity
    Averages(loop,10) = mean(cell_data(:,3));
    Averages(loop,11) = var(cell_data(:,3))/sqrt(length(SD)-1);
    % Direction_cell
    Averages(loop,12) = mean(cell_data(:,4));
    Averages(loop,13) = var(cell_data(:,4))/sqrt(length(SD)-1);
    %DEV
    Averages(loop,14) = mean(summary(:, 8));
    Averages(loop,15) = var(summary(:, 8))/sqrt(length(SD)-1);
    %Signal area
    Averages(loop,16) = mean(mts_area);
    Averages(loop,17) = var(mts_area)/sqrt(length(SD)-1);
    %Aspect ratio
    Averages(loop, 18) = mean(summary(:, 10));
    Averages(loop,19) = var(summary(:, 10))/sqrt(length(SD)-1);
    %Alignment
    Averages(loop,20) = mean(summary(:, 11));
    Averages(loop,21) = var(summary(:, 11))/sqrt(length(SD)-1);
    %Sparseness
    Averages(loop,22) = nansum(summary(:, 12))/length(SD);
    Averages(loop,23) = var(summary(~isnan(summary(:, 12)), 12))/sqrt(length(summary(~isnan(summary(:, 12)), 12))-1);
    close all
end
summary_filename = 'Summary_all.csv';
headers = {'Cell', 'Density', 'sem', 'SD', 'sem', 'Direction_cytoskeleton','sem', 'Area','sem', ...
    'Eccentricity','sem', 'Direction_cell','sem', 'DEV','sem', 'Signal Area','sem',...
    'Aspect ratio','sem', 'Alignment','sem', 'Sperseness','sem'};
cd(sum_dir);
csvwrite_with_headers(summary_filename,Averages,headers);
cd(currdir);

clc
clear variables
close all
