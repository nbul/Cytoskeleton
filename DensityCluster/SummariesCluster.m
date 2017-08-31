%% Writing down summarised data MTSD
cd(MTSD_dir);

summary = zeros(length(SD),8);
for counter2 = 1:length(SD)
    summary(counter2,1) = counter2;
end
summary(:,2) = cell_data(:,3); % cell area
summary(:,3) = cell_data(:,4); % eccentricity
summary(:,4) = cell_data(:,5); % cell orientation
summary(:,5) = SD';
summary(:,6) = mu';
summary(summary(:,6)<0,6) = summary(summary(:,6)<0,6) + 180;
summary(:,7) = 1./sqrt(1-cell_data(:,4).*cell_data(:,4));
summary(:,8) = 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2;
outlier_area = isoutlier(summary(:,2), 'median');
outlier_ecc = outlier_area + isoutlier(summary(:,3), 'median');
outlier_SD = outlier_ecc + isoutlier(summary(:,5), 'median');
outlier_number = length(outlier_SD(outlier_SD ~= 0));
summary(outlier_SD ~= 0,:) = [];
summary_filename = ['Summary_MTSD',num2str(Number),'.csv'];
headers2 = {'Cell', 'Area', 'Eccentricity','Direction_cell', ...
    'SD', 'DEV', 'Elongation', 'Alignment'};
csvwrite_with_headers(summary_filename,summary,headers2);

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
% Sparseness
summary_dens(:,6) = Spars';
% Skewness
summary_dens(:,7) = skew';
% Kurtosis
summary_dens(:,8) = kurt';
% Area
summary_dens(:,9) = cell_data(:,3);
% Eccentricity
summary_dens(:,10) = cell_data(:,4);
summary_dens(isnan(summary_dens(:,3)) == 1,:) = [];
headers3 = {'Cell', 'Signal area', 'Density','Bundling', 'Uniformity', ...
    'Sparseness', 'Skewness', 'Kurtosis', 'Area', 'Eccentricity'};

summary_filename = ['Summary_kmeans_', num2str(Number),'.csv'];
csvwrite_with_headers(summary_filename,summary_dens,headers3);


cd(currdir);