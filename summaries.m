%% Writing down summarised data MTSD


if path==1
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
    summary_filename = [num2str(Number),'_summary_MTSD.csv'];
    headers2 = {'Cell', 'Area', 'Eccentricity','Direction_cell', ...
        'SD', 'DEV', 'Elongation', 'Alignment'};
        cd(MTSD_dir);
    csvwrite_with_headers(summary_filename,summary,headers2);
    cd(currdir);
    Averages(loop,1) = Number;
    % cell area
    Averages(loop,2) = mean(summary(:,2));
    Averages(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
    % eccentricity
    Averages(loop,4) = mean(summary(:,3));
    Averages(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
    % cell orientation
    Averages(loop,6) = mean(summary(:,4));
    Averages(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
    % SD
    Averages(loop,8) = mean(summary(:,5));
    Averages(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
    % Direction cytoskeleton
    Averages(loop,10) = mean(summary(:,6));
    Averages(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
    % cell elongation
    Averages(loop,12) = mean(summary(:, 7));
    Averages(loop,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
    % alignment
    Averages(loop,14) = mean(summary(:, 8));
    Averages(loop,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
    Averages(loop,16) = length(summary(:,1));
    Averages(loop,17) = outlier_number;
else
    summary = zeros(length(mts_area),8);
    for counter2 = 1:length(mts_area)
        summary(counter2,1) = counter2;
    end
    % Signal area
    summary(:,2) = mts_area';
    % Density
    summary(:,3) = mts_density';
    % Bundling
    summary(:,4) = mts_bundling';
    % Uniformity
    summary(:,5) = Uniformity';
    % Sparseness
    summary(:,6) = Spars';
    % Skewness
    summary(:,7) = skew';
    % Kurtosis
    summary(:,8) = kurt';
    % outlier removal
    outlier_2 = isoutlier(summary(:,2), 'median');
    outlier_3 = outlier_2 + isoutlier(summary(:,3), 'median');
    outlier_4 = outlier_3 + isoutlier(summary(:,4), 'median');
    outlier_5 = outlier_4 + isoutlier(summary(:,5), 'median');
    outlier_6 = outlier_5 + isoutlier(summary(:,6), 'median');
    outlier_7 = outlier_6 + isoutlier(summary(:,7), 'median');
    outlier_8 = outlier_7 + isoutlier(summary(:,8), 'median');
    outlier_number = length(outlier_8(outlier_8 ~= 0));
    summary(outlier_8 ~= 0,:) = [];
    summary(isnan(summary(:,3)) == 1,:) = [];
    headers2 = {'Cell', 'Signal area', 'Density','Bundling', 'Uniformity', ...
        'Sparseness', 'Skewness', 'Kurtosis'};
    if method == 1
        cd(edges_dir);
        summary_filename = [num2str(Number),'_summary_edge.csv'];
        csvwrite_with_headers(summary_filename,summary,headers2);
        cd(currdir);
    else
        cd(otsu_dir);
        summary_filename = [num2str(Number),'_summary_otsu.csv'];
        csvwrite_with_headers(summary_filename,summary,headers2);
        cd(currdir);
    end
    Averages(loop,1) = Number;
    % Signal area
    Averages(loop,2) = mean(summary(:,2));
    Averages(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
    % Density
    Averages(loop,4) = mean(summary(:,3));
    Averages(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
    % Bundling
    Averages(loop,6) = mean(summary(:,4));
    Averages(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
    % Uniformity
    Averages(loop,8) = mean(summary(:,5));
    Averages(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
    % Sparseness
    Averages(loop,10) = mean(summary(:,6));
    Averages(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
    % Skewness
    Averages(loop,12) = mean(summary(:, 7));
    Averages(loop,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
    % Kurtosis
    Averages(loop,14) = mean(summary(:, 8));
    Averages(loop,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
    Averages(loop,16) = length(summary(:,1));
    Averages(loop,17) = outlier_number;
end
