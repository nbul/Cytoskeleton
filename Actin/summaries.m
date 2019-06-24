%% Writing down summarised data MTSD

cd(result_dir);

summary = zeros(length(mts_area),3);

summary(:,1) = (1:length(mts_area))';
summary(:,2) = mts_area'; % MT signal area
summary(:,3) = mts_signal'; % MT signal intesity
summary(:,4) = actin_area'; % actin signal area
summary(:,5) = actin_signal'; % actin signal intesity
summary(:,6) = cell_data(:,3); % cell area
summary(:,7) = cell_data(:,4); % eccentricity
summary(:,8) = cell_data(:,5); % cell orientation
summary(:,9) = mu'; % MT direction
summary(:,10) = amu'; % actin direction


summary_filename = ['Summary_MTSD',num2str(loop),'.csv'];

csvwrite_with_headers(summary_filename,summary,headers2);

Pulled_data = [Pulled_data;[(zeros(length(mts_area),1)+1)*loop, summary]];

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
cd(currdir);