clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

sum_dir = [filedir, '/summary'];
dens_dir = [filedir, '/summary/kmeans'];
%% Writing down summarised data MTSD and density embryo by embryo
cd(dens_dir);
files = dir('*.csv');
Averages_dens = zeros(1,1);

for loop=1:numel(files)/2
    
    %% Making summary of density
    
    Density_file = ['Summary_kmeans_', num2str(loop),'.csv'];
    summary_dens = csvread(Density_file,1,0);
    outlier_area = isoutlier(summary_dens(:,10), 'median');
    outlier_ecc = outlier_area + isoutlier(summary_dens(:,11), 'median');
    outlier_number_dens = length(outlier_ecc(outlier_ecc ~= 0));
    summary_dens(outlier_ecc ~= 0,:) = [];
    summary_dens(isnan(summary_dens(:,3)) == 1,:) = [];
    Averages_dens(loop,1) = loop;
    % Signal area
    Averages_dens(loop,2) = mean(summary_dens(:,2));
    Averages_dens(loop,3) = sqrt(var(summary_dens(:,2))/length(summary_dens(:,2)));
    % Density
    Averages_dens(loop,4) = mean(summary_dens(:,3));
    Averages_dens(loop,5) = sqrt(var(summary_dens(:,3))/length(summary_dens(:,3)));
    % Bundling
    Averages_dens(loop,6) = mean(summary_dens(:,4));
    Averages_dens(loop,7) = sqrt(var(summary_dens(:,4))/length(summary_dens(:,4)));
    % Uniformity
    Averages_dens(loop,8) = mean(summary_dens(:,5));
    Averages_dens(loop,9) = sqrt(var(summary_dens(:,5))/length(summary_dens(:,5)));
    % Sparseness
    Averages_dens(loop,10) = mean(summary_dens(:,6));
    Averages_dens(loop,11) = sqrt(var(summary_dens(:,6))/length(summary_dens(:,6)));
    % Skewness
    Averages_dens(loop,12) = mean(summary_dens(:, 7));
    Averages_dens(loop,13) = sqrt(var(summary_dens(:, 7))/length(summary_dens(:, 7)));
    % Kurtosis
    Averages_dens(loop,14) = mean(summary_dens(:, 8));
    Averages_dens(loop,15) = sqrt(var(summary_dens(:, 8))/length(summary_dens(:, 8)));
    
    
    % gaps
    Averages_dens(loop,16) = mean(summary_dens(:,9));
    Averages_dens(loop,17) = sqrt(var(summary_dens(:,9))/length(summary_dens(:,9)));
    % cell number and outliers
    Averages_dens(loop,18) = length(summary_dens(:,1));
    Averages_dens(loop,19) = outlier_number_dens;
    % area
    Averages_dens(loop,20) = mean(summary_dens(:,10));
    Averages_dens(loop,21) = sqrt(var(summary_dens(:,10))/length(summary_dens(:,10)));
    % eccentricity
    Averages_dens(loop,22) = mean(summary_dens(:,11));
    Averages_dens(loop,23) = sqrt(var(summary_dens(:,11))/length(summary_dens(:,11)));
end

cd(sum_dir);

headers_dens = {'Embryo', 'Signal area', 'sem','Density','sem','Bundling','sem',...
    'Uniformity','sem', 'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'Gaps', 'sem', 'Cell number','Outliers', 'Area', 'sem', 'Eccentricity','sem'};
csvwrite_with_headers('Summary_density.csv',Averages_dens,headers_dens);



cd(currdir);
clc
clear variables
close all