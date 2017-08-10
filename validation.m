clear variables
clc

MTnumber = 100;
Eccentricity = 0.9;
shortside = 200;


distribution = 30;
I = 25;
bundling = 1;
method = 1;
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';
%% Setting parameters
parameters = inputdlg({'MT number:','SD:', 'Intensity:','Bundling:','Eccentricity'},...
    'Parameters',1,{num2str(MTnumber), num2str(distribution), num2str(I), num2str(bundling), num2str(Eccentricity)});
% Redefine extension
MTnumber = str2double(parameters{1});
distribution = str2double(parameters{2});
I = str2double(parameters{3});
bundling= str2double(parameters{4});
Eccentricity= str2double(parameters{5});

usedefault = questdlg(strcat('Which threshold methos?'),'Settings','Edge','Otsu','Edge');
if strcmp(usedefault, 'Otsu')
    method = 0;
end

currdir = pwd;
addpath(pwd);
m_added_norm = zeros(45,26);
m_added_norm(:,1) = bincenter;
Value = zeros(45,26);
Value(:,1) = bincenter;
mts_density = zeros(1, 10);
mts_area = zeros(1, 10);
Uniformity = zeros(1, 10);
Spars = zeros(1, 10);
mts_bundling = zeros(1, 10);
kurt = zeros(1, 10);
skew = zeros(1, 10);

long_side = zeros(1, (4*shortside));
for i=(2*shortside):(10*shortside)
    selected_signal = ones(i,shortside);
    to_analyse_o = regionprops(selected_signal, selected_signal,'Eccentricity');
    long_side(i-2*shortside+1) = abs(to_analyse_o.Eccentricity-Eccentricity);
end
[Ecc, Ind] = min(long_side);
longside = Ind + 2*shortside-1;

Length = zeros(1,MTnumber);
result_dir = '/Users/nataliabulgakova/MT-project/Robustness/Densityvalidation';
cd('/Users/nataliabulgakova/MT-project/Robustness/Densityvalidation');
image_dir_name = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_', num2str(method)];
if exist([result_dir,'/', image_dir_name],'dir') == 0
    mkdir(result_dir,image_dir_name);
end
image_dir = [result_dir,'/', image_dir_name];
cd(image_dir);
for k=1:25
    %% Generate random dots within the cell
    rng('shuffle');
    X = randi(shortside, MTnumber,1);
    Y = randi(longside, MTnumber,1);
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
        if Y_temp>0 && Y_temp<=longside
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
        
        X_temp = ceil((longside - b(i))/a(i));
        Y_temp = longside;
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
        if Y_temp>0 && Y_temp<=longside
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
    image = zeros(longside,shortside);
    image_MT_gray = image;
    for i = 1:MTnumber
        image_MT = insertShape(image,'line',intersect(i,:), 'linewidth', 3, 'Color', [I*bundled(i) I*bundled(i) I*bundled(i)]);
        image_MT_gray = image_MT_gray + image_MT(:,:,1);
    end
    
    image_MT_gray = image_MT_gray + 25;
    image_MT_gray(image_MT_gray>255) = 255;
    image_MT_gray = imgaussfilt(image_MT_gray,1);
    
    
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
    background_original = image_MT_gray .* (ones(longside,shortside) - im_bin_c);
    
    %% Density data from signal and background
    selected_signal = ones(longside,shortside);
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
        '_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_', num2str(k),'_', num2str(method), '.tif'];
    print(image, '-dtiff', '-r150', image_filename);
    close all
    
    %% MTSD
     %Apply Sobel Filter over a MTs image to test it
    clear H_full V_full H V M D x y mxd_thr mxd_corrected mxd_indexed
    object_double = image_MT_gray;
    H_full = conv2(object_double,Gx);
    V_full = conv2(object_double,Gy);
    H = H_full(5:longside,5:shortside);
    V = V_full(5:longside,5:shortside);
    M = sqrt(H.^2 + V.^2);
    D = -(180/pi) * atan2(V, H);
    
    [x, y] = size(M);
    p = 1;
    mxd = zeros(1,2);
    
    for j = 2:(y-1)
        for i=2:(x-1)
            %Only directions different to zero are added to table
            if ((M(i,j)) & (M(i+1,j)) & (M(i-1,j)) & (M(i,j+1))  & (M(i,j-1) ~=0) &...
                    (M(i-1,j-1)) & (M(i-1 , j+1)) & (M(i+1 , j-1)) & (M(i+1, j+1))) ~= 0  
                mxd(p,2) = M(i,j); %Second column with magnitudes
                mxd(p,1) = D(i,j); %First column with angles
                p = p + 1;
            end
        end
    end
    max_mxd  = max(mxd(:,2)); %maximum magnitude
    mxd_thr = mxd./repmat([1,max_mxd], length(mxd), 1); %normalised to max magnitude
    
    % Remove all pixels with magnitude less than 22% of maximum
    
    mxd_corrected = mxd_thr(mxd_thr(:,2) >= 0.22,:);
    mxd_corrected(mxd_corrected(:,1) < 0,1) = mxd_corrected(mxd_corrected(:,1) < 0,1) + 180;
    mxd_corrected(:,1) = mxd_corrected(:,1) - 90;
    mxd_corrected(mxd_corrected(:,1) >= 90,1) = 89.9;
    mxd_corrected = sortrows(mxd_corrected,1);
    
    % Make histogram
    [N, bins] = histc(mxd_corrected(:,1),binrange);
    
    mxd_indexed(:,2) = mxd_corrected(:,2);
    mxd_indexed(:,1) = bins(:,1);
    
    
    % Make distribution
    m_added = zeros(45,1);
    for i=1:(length(binrange)-1)
        m_added(i,1) = sum(mxd_indexed(mxd_indexed(:,1) == i,2));
    end
    
    m_added_norm(:,k+1) = m_added/sum(m_added);

    %%Length of microtubules
    for i = 1:MTnumber
        Length(i) = sqrt((intersect(i,1)-intersect(i,3))*(intersect(i,1)-intersect(i,3))...
            + (intersect(i,2)-intersect(i,4))*(intersect(i,2)-intersect(i,4)));
    end
 
    for i=1:(length(binrange)-1)
        Value(i,k+1)=0;
        for l=1:MTnumber
            if angles(l)>=binrange(i) && angles(l)<binrange(i+1)
                Value(i,k+1) = Value(i,k+1) + Length(l)*bundled(l);
            end
        end
    end
    Value(:,k+1) = Value(:,k+1) / sum(Value(:,k+1));
end
cd(currdir);
%% Von Mises
vonmises_fit_dist_sum;
%% Summary
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
% MTSD
summary(:,9) = SD';
% MT direction
mu(mu<0) = 180+mu(mu<0);
summary(:,10) = mu';
m_added_norm = Value;
vonmises_fit_dist_sum;
% MTSD theoretical
summary(:,11) = SD';
% MT theoretical
summary(:,12) = mu'+90;

headers2 = {'Cell', 'Signal area', 'Density','Bundling', 'Uniformity', ...
    'Sparseness', 'Skewness', 'Kurtosis', 'MTSD', 'MT direction', 'MTSD theor', 'MT direction theor'};
summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_',num2str(method),'.csv'];
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
% MTSD
Averages(1,15) = mean(summary(:, 9));
Averages(1,16) = sqrt(var(summary(:, 9))/length(summary(:, 9)));
% MT direction
Averages(1,17) = mean(summary(:, 10));
Averages(1,18) = sqrt(var(summary(:, 10))/length(summary(:, 10)));
% MTSD theoretical
Averages(1,19) = mean(summary(:, 11));
Averages(1,20) = sqrt(var(summary(:, 11))/length(summary(:, 11)));
% MT direction theoretical
Averages(1,21) = mean(summary(:, 12));
Averages(1,22) = sqrt(var(summary(:, 12))/length(summary(:, 12)));
Averages(1,23) = length(summary(:,1));
headers = {'Signal area', 'sem','Density','sem','Bundling','sem', 'Uniformity','sem', ...
    'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'MTSD','sem', 'MT direction','sem', ...
    'MTSD theor','sem', 'MT direction theor','sem', 'Cell number'};
summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_', num2str(method),'_summary.csv'];



cd(result_dir);
csvwrite_with_headers(summary_filename,Averages,headers);
cd(currdir);
clear variables
close all
clc