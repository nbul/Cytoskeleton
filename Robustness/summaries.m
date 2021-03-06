%% Writing down summarised data MTSD

cd(result_dir);

summary = zeros(length(mts_area),3);

summary(:,1) = (1:length(mts_area))';
summary(:,2) = mts_area'; % Signal area
summary(:,3) = mts_signal'; % Signal intesity
summary(:,4) = cell_data(:,3); % cell area
summary(:,5) = cell_data(:,4); % eccentricity
summary(:,6) = cell_data(:,5); % cell orientation
summary(:,7) = SD';
summary(:,8) = mu';
summary(summary(:,8)<0,8) = summary(summary(:,8)<0,8) + 180;
summary(:,9) = 1./sqrt(1-cell_data(:,4).*cell_data(:,4));
summary(:,10) = 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2;
% outlier removal
outlier = zeros(length(mts_area),1);%isoutlier(cell_data(:,3), 'median');
% outlier = outlier + isoutlier(cell_data(:,4), 'median');
% outlier = outlier + isoutlier(summary(:,6), 'median');
outlier_number = length(outlier(outlier ~= 0));
summary(outlier ~= 0,:) = [];

summary_filename = ['Summary_MTSD',num2str(loop),'.csv'];
headers2 = {'Cell', 'Signal area', 'Signal intensity', 'Area', 'Eccentricity','Direction_cell', ...
    'SD', 'DEV', 'Elongation', 'Alignment'};
csvwrite_with_headers(summary_filename,summary,headers2);



%% Averaged values
Averages(loop,1) = loop;
% Signal area
Averages(loop,2) = mean(summary(:,2));
Averages(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Signal intensity
Averages(loop,4) = mean(summary(:,3));
Averages(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
%Cell area
Averages(loop,6) = mean(summary(:,4));
Averages(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
%eccentricity
Averages(loop,8) = mean(summary(:,5));
Averages(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
%Cell orientation
Averages(loop,10) = mean(summary(:,6));
Averages(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
%MTSD
Averages(loop,12) = mean(summary(:,7));
Averages(loop,13) = sqrt(var(summary(:,7))/length(summary(:,7)));
%MT orientation
Averages(loop,14) = mean(summary(:,8));
Averages(loop,15) = sqrt(var(summary(:,8))/length(summary(:,8)));
%elongation
Averages(loop,16) = mean(summary(:,9));
Averages(loop,17) = sqrt(var(summary(:,9))/length(summary(:,9)));
%alignment
Averages(loop,18) = mean(summary(:,10));
Averages(loop,19) = sqrt(var(summary(:,10))/length(summary(:,10)));

Averages(loop,20) = length(summary(:,1));
Averages(loop,21) = outlier_number;
cd(currdir);