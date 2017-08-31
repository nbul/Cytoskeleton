%% Writing down summarised data MTSD

cd(result_dir);
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

summary(isnan(summary(:,3)) == 1,:) = [];
headers2 = {'Cell', 'Signal area', 'Density','Bundling', 'Uniformity', ...
    'Sparseness', 'Skewness', 'Kurtosis'};

summary_filename = [num2str(Number),'_summary_kmeans.csv'];
csvwrite_with_headers(summary_filename,summary,headers2);
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

cd(currdir);