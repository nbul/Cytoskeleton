clear variables
clc

MTnumber = 100;
ratioYX = 2.5;
shortside = 200;
distribution = 30;
I = 25;
bundling = 1;
method = 1;

%% Setting parameters
parameters = inputdlg({'Enter MT number:','Enter SD:', 'Intensity:','Bundling:'},...
    'Parameters',1,{num2str(MTnumber), num2str(distribution), num2str(I), num2str(bundling)});
% Redefine extension
MTnumber = str2double(parameters{1});
distribution = str2double(parameters{2});
I = str2double(parameters{3});
bundling= str2double(parameters{4});

usedefault = questdlg(strcat('Which threshold methos?'),'Settings','Edge','Otsu','Edge');
if strcmp(usedefault, 'Otsu')
    method = 0;
end

currdir = pwd;
addpath(pwd);
mts_density = zeros(1, 10);
mts_area = zeros(1, 10);
Uniformity = zeros(1, 10);
Spars = zeros(1, 10);
mts_bundling = zeros(1, 10);
kurt = zeros(1, 10);
skew = zeros(1, 10);
result_dir = '/Users/nataliabulgakova/MT-project/Robustness/Densityvalidation';
cd('/Users/nataliabulgakova/MT-project/Robustness/Densityvalidation');
image_dir_name = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_', num2str(method)];
if exist([result_dir,'/', image_dir_name],'dir') == 0
    mkdir(result_dir,image_dir_name);
end
image_dir = [result_dir,'/', image_dir_name];
cd(image_dir);
for k=1:25
    %% Generate random dots within the cell
    rng('shuffle');
    X = randi(shortside, MTnumber,1);
    Y = randi(ratioYX*shortside, MTnumber,1);
    %% Generate angles with a given distribution
    angles = normrnd(0,distribution, [MTnumber 1]);
    for i=1:MTnumber
        if angles(i, 1) <= -90
            angles(i, 1) = angles(i, 1) + 180;
        elseif angles(i, 1) > 90
            angles(i, 1) = angles(i, 1) - 180;
        end
    end
    
    bundled = randi(bundling, MTnumber,1);
    %% Line parameters and start/end points
    a = 1./tand(angles);
    b = Y - a.*X;
    intersect=zeros(MTnumber,4);
    l=zeros(MTnumber,1);
    for i=1:MTnumber
        l(i)=0;
        X_temp = 1;
        Y_temp = ceil(b(i)+a(i));
        if Y_temp>0 && Y_temp<=ratioYX*shortside
            l(i)=1;
            intersect(i,1) = X_temp;
            intersect(i,2) = Y_temp;
        end
        
        X_temp = ceil(1 - b(i)/a(i));
        Y_temp = 1;
        if X_temp>0 && X_temp<=shortside
            if l(i) == 0
                l(i) = 1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            else
                l(i) = l(i)+1;
                intersect(i,3) = X_temp;
                intersect(i,4) = Y_temp;
            end
        end
        
        X_temp = ceil((ratioYX*shortside - b(i))/a(i));
        Y_temp = ratioYX*shortside;
        if X_temp>0 && X_temp<=shortside
            if l(i) == 0
                l(i) = 1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            else
                l(i) = l(i)+1;
                intersect(i,3) = X_temp;
                intersect(i,4) = Y_temp;
            end
        end
        
        X_temp = shortside;
        Y_temp = a(i)*shortside +b(i);
        if Y_temp>0 && Y_temp<=ratioYX*shortside
            if l(i) == 0
                l(i) = 1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            else
                l(i) = l(i)+1;
                intersect(i,3) = X_temp;
                intersect(i,4) = Y_temp;
            end
        end
    end
    
    %% Draw lines
    image = zeros(ratioYX*shortside,shortside);
    image_MT_gray = image;
    for i = 1:MTnumber
        image_MT = insertShape(image,'line',intersect(i,:), 'linewidth', 5, 'Color', [I*bundled(i) I*bundled(i) I*bundled(i)]);
        image_MT_gray = image_MT_gray + image_MT(:,:,1);
    end
    
    image_MT_gray = image_MT_gray + 25;
    image_MT_gray(image_MT_gray>255) = 255;
    image_MT_gray = imgaussfilt(image_MT_gray,2);
    
    
    image_edges = edge(image_MT_gray, 'Canny');
    signal_edges = image_MT_gray .* image_edges;
    %% Edge detection
    if method == 1
        [A1, A2] = histcounts(signal_edges(signal_edges>0));
        bins = (min(A2)+((A2(2)-A2(1))/2)):((A2(2)-A2(1))):(max(A2)-((A2(2)-A2(1))/2));
        [XOut,YOut] = prepareCurveData(bins,A1);
        fo = fitoptions('gauss1', 'Lower', [0 min(A2) 0], 'Upper', [Inf max(A2) Inf]);
        [threshold, gof_edges] = fit(XOut, YOut, 'gauss1', fo);
        im_bin_c = imbinarize(image_MT_gray,threshold.b1*0.7);
    else
        im_bin_c = imbinarize(imadjust(image_MT_gray/255),graythresh(imadjust(image_MT_gray/255))*0.7);
    end
    %% Generate Cell Masks.
    signal_original = image_MT_gray .* im_bin_c;
    background_original = image_MT_gray .* (ones(ratioYX*shortside,shortside) - im_bin_c);
    
    %% Density data from signal and background
    selected_signal = ones(ratioYX*shortside,shortside);
    to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
    to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');
    
    % Relative to background and signal area
    sum_pixvalues_o = sum(to_analyse_o.PixelValues(:,1));
    sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
    num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
    num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
    
    mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    
    % Signal Area
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
    
    % Relative to edges intensity
    if max(to_analyse_c.PixelValues)~= 0
        mts_bundling(k) = (mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))-...
            (0.7*mean(signal_edges(signal_edges>0)))) / ...
            (0.7*mean(signal_edges(signal_edges>0)));
        %         mts_bundling(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        %             (sum_pixvalues_back_o / num_pixvalues_back_c));
        %         mts_bundling(k) = mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))...
        %             /min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
    else
        mts_bundling(k) = 0;
    end
    
    % Uniformity
    Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
        (to_analyse_all.PixelValues + mean(to_analyse_all.PixelValues)))/length(to_analyse_all.PixelValues));
    %Sparseness
    Spars(k) = calcSparseness(to_analyse_o.PixelValues/mean(to_analyse_o.PixelValues(to_analyse_o.PixelValues>0)),1);
    
    % Kurtosis and skewness
    signal = to_analyse_all.PixelValues - (0.7*mean(signal_edges(signal_edges>0)));
    kurt(k) = kurtosis(signal);
    skew(k) = skewness(signal);
    image = figure;
    imshow(image_MT_gray, [0 255]);
    image_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',num2str(I),...
        '_bund',num2str(bundling),'_', num2str(k),'_', num2str(method), '.tif'];
    print(image, '-dtiff', '-r150', image_filename);
    close all
end
cd(currdir);
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
summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_',num2str(method),'.csv'];
cd(result_dir);
csvwrite_with_headers(summary_filename,summary,headers2);
cd(currdir);
% Signal area
Averages(1,1) = mean(summary(:,2));
Averages(1,2) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Density
Averages(1,3) = mean(summary(:,3));
Averages(1,4) = sqrt(var(summary(:,3))/length(summary(:,3)));
% Bundling
Averages(1,5) = mean(summary(:,4));
Averages(1,6) = sqrt(var(summary(:,4))/length(summary(:,4)));
% Uniformity
Averages(1,7) = mean(summary(:,5));
Averages(1,8) = sqrt(var(summary(:,5))/length(summary(:,5)));
% Sparseness
Averages(1,9) = mean(summary(:,6));
Averages(1,10) = sqrt(var(summary(:,6))/length(summary(:,6)));
% Skewness
Averages(1,11) = mean(summary(:, 7));
Averages(1,12) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
% Kurtosis
Averages(1,13) = mean(summary(:, 8));
Averages(1,14) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
Averages(1,15) = length(summary(:,1));
headers = {'Signal area', 'sem','Density','sem','Bundling','sem', 'Uniformity','sem', ...
    'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem','Cell number'};
summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_', num2str(method),'_summary.csv'];
cd(result_dir);
csvwrite_with_headers(summary_filename,Averages,headers);
cd(currdir);
clear variables
close all
clc