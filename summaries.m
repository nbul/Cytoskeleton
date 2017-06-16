%% Writing down summarised data

summary = zeros(length(SD),12);
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

summary_filename = [num2str(Number),'_summary.csv'];
headers = {'Cell', 'Density', 'SD', 'Direction_cytoskeleton','Area', 'Eccentricity',...
    'Dorection_cell','DEV', 'Signal Area','Aspect ratio','Alignment', 'Sparseness','Bundling'};
cd(sum_dir);
csvwrite_with_headers(summary_filename,summary,headers);

%% writing down averaged by embryo data
Averages(loop,1) = Number;
% Density
Averages(loop,2) = nanmean(mts_density);
Averages(loop,3) = sqrt(var(mts_density(~isnan(mts_density)))/length(mts_density(~isnan(mts_density))));
%SD
Averages(loop,4) = mean(SD);
Averages(loop,5) = sqrt(var(SD)/length(SD));
%Direction cytoskeleton
Averages(loop,6) = mean(mu);
Averages(loop,7) = sqrt(var(mu)/length(SD));
%Area
Averages(loop,8) = mean(cell_data(:,3));
Averages(loop,9) = sqrt(var(cell_data(:,3))/length(SD));
% Eccentricity
Averages(loop,10) = mean(cell_data(:,4));
Averages(loop,11) = sqrt(var(cell_data(:,4))/length(SD));
% Direction_cell
Averages(loop,12) = mean(cell_data(:,5));
Averages(loop,13) = sqrt(var(cell_data(:,5))/length(SD));
%DEV
Averages(loop,14) = mean(summary(:, 8));
Averages(loop,15) = sqrt(var(summary(:, 8))/length(SD));
%Signal area
Averages(loop,16) = mean(mts_area);
Averages(loop,17) = sqrt(var(mts_area)/length(SD));
%Aspect ratio
Averages(loop,18) = mean(summary(:, 10));
Averages(loop,19) = sqrt(var(summary(:, 10))/length(SD));
%Alignment
Averages(loop,20) = mean(summary(:, 11));
Averages(loop,21) = sqrt(var(summary(:, 11))/length(SD));
%Sparseness
Averages(loop,22) = nansum(summary(:, 12))/length(SD);
Averages(loop,23) = sqrt(var(summary(~isnan(summary(:, 12)), 12))/...
    length(summary(~isnan(summary(:, 12)), 12)));
%Bundling
%SD
Averages(loop,24) = mean(summary(:,13));
Averages(loop,25) = sqrt(var(summary(:,13))/length(SD));

%number of cells
Averages(loop,26) = length(SD);
