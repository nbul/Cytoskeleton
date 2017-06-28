%% Writing down summarised data

summary = zeros(length(SD),14);
for counter2 = 1:length(SD)
    summary(counter2,1) = counter2;
end
summary(:,2) = mts_density';
summary(:,3) = SD';
summary(:,4) = mu';
summary(:,5) = cell_data(:,3);
summary(:,6) = cell_data(:,4);
summary(:,7) = cell_data(:,5);
for counter2 = 1:length(SD)
    if (abs(summary(counter2, 4)-summary(counter2, 7)) >= 90)
        summary(counter2, 8) = 180 - abs(summary(counter2, 4)-summary(counter2, 7));
    else
        summary(counter2, 8) = abs(summary(counter2, 4)-summary(counter2, 7));
    end
end
summary(:,9) = mts_area';
summary(:, 10) = 1./sqrt(1-cell_data(:,4).*cell_data(:,4));
summary(:,11) = 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2;
summary(:,12) = Spars';
summary(:,13) = mts_bundling';
summary(:,14) = Uniformity';

summary(summary(:,3)>=90,:) = [];
summary_filename = [num2str(Number),'_summary.csv'];
headers = {'Cell', 'Density', 'SD', 'Direction_cytoskeleton','Area', 'Eccentricity',...
    'Dorection_cell','DEV', 'Signal Area','Aspect ratio','Alignment', 'Sparseness',...
    'Bundling','Uniformity'};
cd(sum_dir);
csvwrite_with_headers(summary_filename,summary,headers);

%% writing down averaged by embryo data
Averages(loop,1) = Number;
% Density
Averages(loop,2) = nanmean(summary(:,2));
Averages(loop,3) = sqrt(var(summary(~isnan(summary(:, 2)), 2))/...
    length(mts_density(~isnan(summary(:,2)))));
%SD
Averages(loop,4) = mean(summary(:,3));
Averages(loop,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
%Direction cytoskeleton
Averages(loop,6) = mean(summary(:,4));
Averages(loop,7) = sqrt(var(summary(:,4))/length(summary(:,4)));
%Area
Averages(loop,8) = mean(summary(:,5));
Averages(loop,9) = sqrt(var(summary(:,5))/length(summary(:,5)));
% Eccentricity
Averages(loop,10) = mean(summary(:,6));
Averages(loop,11) = sqrt(var(summary(:,6))/length(summary(:,6)));
% Direction_cell
Averages(loop,12) = mean(summary(:,7));
Averages(loop,13) = sqrt(var(summary(:,7))/length(summary(:,7)));
%DEV
Averages(loop,14) = mean(summary(:, 8));
Averages(loop,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
%Signal area
Averages(loop,16) = mean(summary(:,9));
Averages(loop,17) = sqrt(var(summary(:,9))/length(summary(:,9)));
%Aspect ratio
Averages(loop,18) = mean(summary(:, 10));
Averages(loop,19) = sqrt(var(summary(:, 10))/length(summary(:,10)));
%Alignment
Averages(loop,20) = mean(summary(:, 11));
Averages(loop,21) = sqrt(var(summary(:, 11))/length(summary(:, 11)));
%Sparseness
Averages(loop,22) = nansum(summary(:, 12))/length(SD);
Averages(loop,23) = sqrt(var(summary(~isnan(summary(:, 12)), 12))/...
    length(summary(~isnan(summary(:, 12)), 12)));
%Bundling
Averages(loop,24) = mean(summary(:,13));
Averages(loop,25) = sqrt(var(summary(:,13))/length(summary(:,13)));

%Uniformity
Averages(loop,26) = mean(summary(:,14));
Averages(loop,27) = sqrt(var(summary(:,14))/length(summary(:,14)));

%number of cells
Averages(loop,28) = length(summary(:,13));
