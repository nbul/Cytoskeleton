%% Writing down summarised data MTSD

cd(result_dir);

summary = zeros(length(mts_area),3);

summary(:,1) = (1:length(mts_area))';
summary(:,2) = mts_area'; % Signal area
summary(:,3) = cell_data(:,3); % cell area
summary(:,4) = cell_data(:,4); % eccentricity
summary(:,5) = cell_data(:,5); % cell orientation
summary(:,6) = SD';
summary(:,7) = mu';
summary(summary(:,7)<0,7) = summary(summary(:,7)<0,7) + 180;
summary(:,8) = 1./sqrt(1-cell_data(:,4).*cell_data(:,4));
summary(:,9) = 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2;
% outlier removal
outlier = zeros(length(mts_area),1);%isoutlier(cell_data(:,3), 'median');
% outlier = outlier + isoutlier(cell_data(:,4), 'median');
% outlier = outlier + isoutlier(summary(:,6), 'median');
outlier_number = length(outlier(outlier ~= 0));
summary(outlier ~= 0,:) = [];

summary_filename = ['Summary_MTSD',num2str(loop),'.csv'];
headers2 = {'Cell', 'Signal', 'Area', 'Eccentricity','Direction_cell', ...
    'SD', 'DEV', 'Elongation', 'Alignment'};
csvwrite_with_headers(summary_filename,summary,headers2);



%% Averaged values
Averages(loop,1) = loop;
% Signal area
Averages(loop,2) = mean(summary(:,2));
Averages(loop,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Cell area
Averages(loop,4) = mean(summary(:,3));
Averages(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
%eccentricity
Averages(loop,6) = mean(summary(:,4));
Averages(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
%Cell orientation
Averages(loop,8) = mean(summary(:,5));
Averages(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
%MTSD
Averages(loop,10) = mean(summary(:,6));
Averages(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
%MT orientation
Averages(loop,12) = mean(summary(:,7));
Averages(loop,13) = sqrt(var(summary(:,7))/length(summary(:,7)));
%elongation
Averages(loop,14) = mean(summary(:,8));
Averages(loop,15) = sqrt(var(summary(:,8))/length(summary(:,8)));
%alignment
Averages(loop,16) = mean(summary(:,9));
Averages(loop,17) = sqrt(var(summary(:,9))/length(summary(:,9)));

Averages(loop,18) = length(summary(:,1));
Averages(loop,19) = outlier_number;
cd(currdir);