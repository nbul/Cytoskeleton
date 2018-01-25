%% Writing down summarised data MTSD

cd(result_dir);

summary = zeros(length(mts_area),3);

summary(:,1) = (1:length(mts_area))';

% Signal area
summary(:,2) = mts_area';
% Bundling
summary(:,3) = mts_bundling';

% outlier removal
outlier = isoutlier(cell_data(:,3), 'median');
outlier = outlier + isoutlier(cell_data(:,4), 'median');
outlier_number = length(outlier(outlier ~= 0));
summary(outlier ~= 0,:) = [];

headers2 = {'Cell', 'Signal area', 'Bundling'};
summary_filename = [num2str(loop),'_summary.csv'];
csvwrite_with_headers(summary_filename,summary,headers2);
%% Averaged values
Averages(loop,1) = loop;
% Signal area
Averages(loop,2) = mean(summary(:,2));
Averages(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Bundling
Averages(loop,4) = mean(summary(:,3));
Averages(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
Averages(loop,6) = length(summary(:,1));
Averages(loop,7) = outlier_number;
cd(currdir);