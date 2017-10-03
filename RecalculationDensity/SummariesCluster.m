%% Writing down summarised data Density
cd(dens_dir);
summary_dens = zeros(length(mts_area),8);
for counter2 = 1:length(mts_area)
    summary_dens(counter2,1) = counter2;
end
% Signal area
summary_dens(:,2) = mts_area';
% Density
summary_dens(:,3) = mts_density';
% Bundling
summary_dens(:,4) = mts_bundling';
% Uniformity
summary_dens(:,5) = Uniformity';
% UNAAD
summary_dens(:,6) = UNAAD';
% Sparseness
summary_dens(:,7) = Spars';
% Skewness
summary_dens(:,8) = skew';
% Kurtosis
summary_dens(:,9) = kurt';
% Gaps
summary_dens(:,10) = space';
% Area
summary_dens(:,11) = cell_data(:,3);
% Eccentricity
summary_dens(:,12) = cell_data(:,4);
summary_dens(isnan(summary_dens(:,3)) == 1,:) = [];
headers3 = {'Cell', 'Signal area', 'Density','Bundling', 'Uniformity', 'UNAAD', ...
    'Sparseness', 'Skewness', 'Kurtosis', 'Gaps', 'Area', 'Eccentricity'};

summary_filename = ['Summary_kmeans_', num2str(Number),'_post.csv'];
csvwrite_with_headers(summary_filename,summary_dens,headers3);

outlier_area = isoutlier(summary_dens(:,11), 'median');
outlier_ecc = outlier_area + isoutlier(summary_dens(:,12), 'median');
outlier_number_dens = length(outlier_ecc(outlier_ecc ~= 0));
summary_dens(outlier_ecc ~= 0,:) = [];
summary_dens(isnan(summary_dens(:,3)) == 1,:) = [];
Averages_dens(loop,1) = Number;
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
% UNAAD
Averages_dens(loop,10) = mean(summary_dens(:,6));
Averages_dens(loop,11) = sqrt(var(summary_dens(:,6))/length(summary_dens(:,6)));
% Sparseness
Averages_dens(loop,12) = mean(summary_dens(:, 7));
Averages_dens(loop,13) = sqrt(var(summary_dens(:, 7))/length(summary_dens(:, 7)));
% Skewness
Averages_dens(loop,14) = mean(summary_dens(:, 8));
Averages_dens(loop,15) = sqrt(var(summary_dens(:, 8))/length(summary_dens(:, 8)));
% Kurtosis
Averages_dens(loop,16) = mean(summary_dens(:,9));
Averages_dens(loop,17) = sqrt(var(summary_dens(:,9))/length(summary_dens(:,9)));
% gaps
Averages_dens(loop,18) = length(summary_dens(:,1));
Averages_dens(loop,19) = outlier_number_dens;
% cell number and outliers
Averages_dens(loop,20) = mean(summary_dens(:,10));
Averages_dens(loop,21) = sqrt(var(summary_dens(:,10))/length(summary_dens(:,10)));
% area
Averages_dens(loop,22) = mean(summary_dens(:,11));
Averages_dens(loop,23) = sqrt(var(summary_dens(:,11))/length(summary_dens(:,11)));
% eccentricity
Averages_dens(loop,24) = mean(summary_dens(:,11));
Averages_dens(loop,25) = sqrt(var(summary_dens(:,11))/length(summary_dens(:,11)));

cd(currdir);