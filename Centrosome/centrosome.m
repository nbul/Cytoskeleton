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
cd(b_dir_control);
files_control = dir('*.tif');
cd(currdir);

%% Determining threshold
threshold = zeros(length(files_control),1);
for loop=1:length(files_control)
    cd(cent_dir_control);
    Cs_file = [num2str(loop),'.tif'];
    Image = imread(Cs_file);
    Image = imgaussfilt(Image,2);
    image_original_double = im2double(Image);
    threshold(loop) = graythresh(imadjust(image_original_double))*2;
%     im_bin_c = imbinarize(imadjust(image_original_double),0.8566);
%     imshow(bwareaopen(im_bin_c,10));
end
thr = mean(threshold);
%% Analysis

headers2 = {'Cell', 'Original number', 'Area', 'Eccentricity','Direction_cell', ...
    'Number spots', 'Total intensity', 'Distance'};

headers1 = {'Embryo', 'Original number', 'sem', 'Area', 'sem', 'Eccentricity', 'sem','Direction_cell', 'sem', ...
    'Number spots', 'sem', 'Total intensity', 'sem', 'Distance', 'sem','Number cells'};

control;
mutant;

cd(currdir);
clc
clear variables
close all