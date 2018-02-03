%% Writing down summarised data MTSD
cd(result_dir_MTSD)
summary_MTSD = zeros(length(SD),8);
for counter2 = 1:length(SD)
    summary_MTSD(counter2,1) = counter2;
end
summary_MTSD(:,2) = cell_data(:,3); % cell area
summary_MTSD(:,3) = cell_data(:,4); % eccentricity
summary_MTSD(:,4) = cell_data(:,5); % cell orientation
summary_MTSD(:,5) = SD';
summary_MTSD(:,6) = mu';
summary_MTSD(summary_MTSD(:,6)<0,6) = summary_MTSD(summary_MTSD(:,6)<0,6) + 180;
summary_MTSD(:,7) = 1./sqrt(1-cell_data(:,4).*cell_data(:,4));
summary_MTSD(:,8) = 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2;
outlier_area = isoutlier(summary_MTSD(:,2), 'median');
outlier_ecc = outlier_area + isoutlier(summary_MTSD(:,3), 'median');
outlier_SD = outlier_ecc + isoutlier(summary_MTSD(:,5), 'median');
outlier_number = length(outlier_SD(outlier_SD ~= 0));
summary_MTSD(outlier_SD ~= 0,:) = [];
summary_MTSD_filename = [num2str(Number),'_summary_MTSD.csv'];
headers2_MTSD = {'Cell', 'Area', 'Eccentricity','Direction_cell', ...
    'SD', 'DEV', 'Elongation', 'Alignment'};
csvwrite_with_headers(summary_MTSD_filename,summary_MTSD,headers2_MTSD);
%% Writing down averaged data MTSD
Averages_MTSD(loop,1) = Number;
% cell area
Averages_MTSD(loop,2) = mean(summary_MTSD(:,2));
Averages_MTSD(loop,3) = sqrt(var(summary_MTSD(:,2))/length(summary_MTSD(:,2)));
% eccentricity
Averages_MTSD(loop,4) = mean(summary_MTSD(:,3));
Averages_MTSD(loop,5) = sqrt(var(summary_MTSD(:,3))/length(summary_MTSD(:,3)));
% cell orientation
Averages_MTSD(loop,6) = mean(summary_MTSD(:,4));
Averages_MTSD(loop,7) = sqrt(var(summary_MTSD(:,4))/length(summary_MTSD(:,4)));
% SD
Averages_MTSD(loop,8) = mean(summary_MTSD(:,5));
Averages_MTSD(loop,9) = sqrt(var(summary_MTSD(:,5))/length(summary_MTSD(:,5)));
% Direction cytoskeleton
Averages_MTSD(loop,10) = mean(summary_MTSD(:,6));
Averages_MTSD(loop,11) = sqrt(var(summary_MTSD(:,6))/length(summary_MTSD(:,6)));
% cell elongation
Averages_MTSD(loop,12) = mean(summary_MTSD(:, 7));
Averages_MTSD(loop,13) = sqrt(var(summary_MTSD(:, 7))/length(summary_MTSD(:, 7)));
% alignment
Averages_MTSD(loop,14) = mean(summary_MTSD(:, 8));
Averages_MTSD(loop,15) = sqrt(var(summary_MTSD(:, 8))/length(summary_MTSD(:, 8)));
Averages_MTSD(loop,16) = length(summary_MTSD(:,1));
Averages_MTSD(loop,17) = outlier_number;


%% Summarised data for density
cd(result_dir_Dens);
summary_Dens = zeros(length(Uniformity),7);
for counter2 = 1:length(Uniformity)
    summary_Dens(counter2,1) = counter2;
end
% Uniformity
summary_Dens(:,2) = Uniformity';
% Sparseness
summary_Dens(:,3) = Spars';
% Skewness
summary_Dens(:,4) = skew';
% Kurtosis
summary_Dens(:,5) = kurt';
% Sdr
summary_Dens(:,6) = Sdr';
% Sdq
summary_Dens(:,7) = Sdq';

% outlier removal
outlier_2 = isoutlier(summary_Dens(:,7), 'median');
%     outlier_3 = outlier_2 + isoutlier(summary(:,3), 'median');
%     outlier_4 = outlier_3 + isoutlier(summary(:,4), 'median');
%     outlier_5 = outlier_4 + isoutlier(summary(:,5), 'median');
%     outlier_6 = outlier_5 + isoutlier(summary(:,6), 'median');
%     outlier_7 = outlier_6 + isoutlier(summary(:,7), 'median');
%     outlier_8 = outlier_7 + isoutlier(summary(:,8), 'median');
outlier_8 = outlier_2 + isoutlier(cell_data(:,3), 'median');
outlier_number = length(outlier_8(outlier_8 ~= 0));
summary_Dens(outlier_8 ~= 0,:) = [];
summary_Dens(isnan(summary_Dens(:,3)) == 1,:) = [];
headers2_Dens = {'Cell', 'Uniformity', 'Sparseness', 'Skewness', 'Kurtosis',...
    'Sdr', 'Sdq'};

summary_filename_Dens = [num2str(Number),'_summary_Dens.csv'];
csvwrite_with_headers(summary_filename_Dens,summary_Dens,headers2_Dens);

%% Writing averaged data for density
Averages_Dens(loop,1) = Number;
% Uniformity
Averages_Dens(loop,2) = mean(summary_Dens(:,2));
Averages_Dens(loop,3) = sqrt(var(summary_Dens(:,2))/length(summary_Dens(:,2)));
% Sparseness
Averages_Dens(loop,4) = mean(summary_Dens(:,3));
Averages_Dens(loop,5) = sqrt(var(summary_Dens(:,3))/length(summary_Dens(:,3)));
% Skewness
Averages_Dens(loop,6) = mean(summary_Dens(:,4));
Averages_Dens(loop,7) = sqrt(var(summary_Dens(:,4))/length(summary_Dens(:,4)));
% Kurtosis
Averages_Dens(loop,8) = mean(summary_Dens(:,5));
Averages_Dens(loop,9) = sqrt(var(summary_Dens(:,5))/length(summary_Dens(:,5)));
% Sdr
Averages_Dens(loop,10) = mean(summary_Dens(:,6));
Averages_Dens(loop,11) = sqrt(var(summary_Dens(:,6))/length(summary_Dens(:,6)));
% Sdq
Averages_Dens(loop,12) = mean(summary_Dens(:, 7));
Averages_Dens(loop,13) = sqrt(var(summary_Dens(:, 7))/length(summary_Dens(:, 7)));

Averages_Dens(loop,14) = length(summary_Dens(:,1));
Averages_Dens(loop,15) = outlier_number;
cd(currdir);