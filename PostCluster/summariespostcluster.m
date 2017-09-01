clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
image_dir =[filedir, '/cytoskeleton_average'];
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];
b_dir = [filedir, '/borders'];
sum_dir = [filedir, '/summary'];
dens_dir = [filedir, '/summary/kmeans'];
MTSD_dir = [filedir, '/summary/MTSD'];
%% Writing down summarised data MTSD
cd(dens_dir);
files = dir('*.csv');
Averages_dens = zeros(1,1);
Averages_MTSD = zeros(1,1);
for loop=1:numel(files)
    %% reading files
    cd(image_dir);
    clear Name Number Actin_file Image_actin Path  Image_borders
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
    
    %% Open MTs image, adjust and generate cell masks
    imageprocessing;
    %% Making summary of density
    cd(MTSD_dir);
    MTSD_file = ['Summary_MTSD',num2str(Number),'.csv'];
    summary_dens = csvread(MTSD_file,1,0);
    outlier_area = isoutlier(summary_dens(:,2), 'median');
    outlier_ecc = outlier_area + isoutlier(summary_dens(:,3), 'median');
    outlier_SD = outlier_ecc + isoutlier(summary_dens(:,5), 'median');
    outlier_number = length(outlier_SD(outlier_SD ~= 0));
    summary_dens(outlier_SD ~= 0,:) = [];
    Averages_MTSD(loop,1) = Number;
    % cell area
    Averages_MTSD(loop,2) = mean(summary(:,2));
    Averages_MTSD(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
    % eccentricity
    Averages_MTSD(loop,4) = mean(summary(:,3));
    Averages_MTSD(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
    % cell orientation
    Averages_MTSD(loop,6) = mean(summary(:,4));
    Averages_MTSD(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
    % SD
    Averages_MTSD(loop,8) = mean(summary(:,5));
    Averages_MTSD(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
    % Direction cytoskeleton
    Averages_MTSD(loop,10) = mean(summary(:,6));
    Averages_MTSD(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
    % cell elongation
    Averages_MTSD(loop,12) = mean(summary(:, 7));
    Averages_MTSD(loop,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
    % alignment
    Averages_MTSD(loop,14) = mean(summary(:, 8));
    Averages_MTSD(loop,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
    Averages_MTSD(loop,16) = length(summary(:,1));
    Averages_MTSD(loop,17) = outlier_number;
    
    %% averages density
    cd(dens_dir);
    Density_file = ['Summary_kmeans_', num2str(Number),'.csv'];
    summary_MTSD = csvread(MTSD_file,1,0);
    outlier_area = isoutlier(summary_MTSD(:,9), 'median');
    outlier_ecc = outlier_area + isoutlier(summary_MTSD(:,10), 'median');
    outlier_number = length(outlier_ecc(outlier_ecc ~= 0));
    summary_MTSD(outlier_ecc ~= 0,:) = [];
    summary_MTSD(isnan(summary(:,3)) == 1,:) = [];
    Averages_dens(loop,1) = Number;
    % Signal area
    Averages_dens(loop,2) = mean(summary(:,2));
    Averages_dens(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
    % Density
    Averages_dens(loop,4) = mean(summary(:,3));
    Averages_dens(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
    % Bundling
    Averages_dens(loop,6) = mean(summary(:,4));
    Averages_dens(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
    % Uniformity
    Averages_dens(loop,8) = mean(summary(:,5));
    Averages_dens(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
    % Sparseness
    Averages_dens(loop,10) = mean(summary(:,6));
    Averages_dens(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
    % Skewness
    Averages_dens(loop,12) = mean(summary(:, 7));
    Averages_dens(loop,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
    % Kurtosis
    Averages_dens(loop,14) = mean(summary(:, 8));
    Averages_dens(loop,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
    Averages_dens(loop,16) = length(summary(:,1));
    Averages_dens(loop,17) = outlier_number;
    % area
    Averages_dens(loop,18) = mean(summary(:,9));
    Averages_dens(loop,19) = sqrt(var(summary(:,9))/length(summary(:,9)));
    % eccentricity
    Averages_dens(loop,20) = mean(summary(:,10));
    Averages_dens(loop,21) = sqrt(var(summary(:,10))/length(summary(:,10)));
end

cd(sum_dir);
headers_MTSD = {'Embryo', 'Area','sem',  'Eccentricity','sem', 'Direction_cell','sem',...
    'SD', 'sem', 'Direction_cytoskeleton','sem', 'Aspect ratio','sem', 'Alignment','sem',...
    'Cell number','Outliers'};
csvwrite_with_headers('Summary_MTSD.csv',Averages_MTSD,headers_MTSD);

headers_dens = {'Embryo', 'Signal area', 'sem','Density','sem','Bundling','sem',...
    'Uniformity','sem', 'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'Cell number','Outliers', 'Area', 'sem', 'Eccentricity','sem'};
csvwrite_with_headers('Summary_density.csv',Averages_dens,headers_dens);

cd(currdir);
clc
clear variables
close all