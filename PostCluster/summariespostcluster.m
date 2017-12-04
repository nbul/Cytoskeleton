clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
image_dir =[filedir, '/cytoskeleton_average'];
image_dir2 =[filedir, '/cytoskeleton'];
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];
b_dir = [filedir, '/borders'];
dist_dir = [filedir, '/distribution'];
sum_dir = [filedir, '/summary'];
dens_dir = [filedir, '/summary/kmeans'];
MTSD_dir = [filedir, '/summary/MTSD'];
%% Writing down summarised data MTSD and density embryo by embryo
cd(MTSD_dir);
files = dir('*.csv');
Averages_dens = zeros(1,1);
Averages_MTSD = zeros(1,1);
for loop=1:numel(files)
    
    %% reading files
    cd(image_dir);
    clear Name Number Actin_file Image_actin Path  Image_borders
    Number = loop;
    Actin_file = [num2str(Number),'.tif'];
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('tracked_bd.png');
    [im_x, im_y] = size(Image2);
    
    %% Collect data about cells and boundaries
    borders;
    
    %% Open MTs image, adjust and generate cell masks
    imageprocessing;
    %% Making summary of density
    cd(MTSD_dir);
    MTSD_file = ['Summary_MTSD',num2str(Number),'.csv'];
    summary_MTSD = csvread(MTSD_file,1,0);
    outlier_area = isoutlier(summary_MTSD(:,2), 'median');
    outlier_ecc = outlier_area + isoutlier(summary_MTSD(:,3), 'median');
    outlier_SD = outlier_ecc + isoutlier(summary_MTSD(:,5), 'median');
    outlier_number_MTSD = length(outlier_SD(outlier_SD ~= 0));
    summary_MTSD(outlier_SD ~= 0,:) = [];
    Averages_MTSD(loop,1) = Number;
    % cell area
    Averages_MTSD(loop,2) = mean(summary_MTSD(:,2));
    Averages_MTSD(loop,3) = sqrt(var(summary_MTSD(:,2))/length(summary_MTSD(:,2)));
    % eccentricity
    Averages_MTSD(loop,4) = mean(summary_MTSD(:,3));
    Averages_MTSD(loop,5) = sqrt(var(summary_MTSD(:,3))/length(summary_MTSD(:,3)));
    % cell orientation
    Averages_MTSD(loop,6) = mean(summary_MTSD(:,4));
    Averages_MTSD(loop,7) = sqrt(var(summary_MTSD(:,4))/length(summary_MTSD(:,4)));
    % SD
    Averages_MTSD(loop,8) = mean(summary_MTSD(:,5));
    Averages_MTSD(loop,9) = sqrt(var(summary_MTSD(:,5))/length(summary_MTSD(:,5)));
    % Direction cytoskeleton
    Averages_MTSD(loop,10) = mean(summary_MTSD(:,6));
    Averages_MTSD(loop,11) = sqrt(var(summary_MTSD(:,6))/length(summary_MTSD(:,6)));
    % cell elongation
    Averages_MTSD(loop,12) = mean(summary_MTSD(:, 7));
    Averages_MTSD(loop,13) = sqrt(var(summary_MTSD(:, 7))/length(summary_MTSD(:, 7)));
    % alignment
    Averages_MTSD(loop,14) = mean(summary_MTSD(:, 8));
    Averages_MTSD(loop,15) = sqrt(var(summary_MTSD(:, 8))/length(summary_MTSD(:, 8)));
    Averages_MTSD(loop,16) = length(summary_MTSD(:,1));
    Averages_MTSD(loop,17) = outlier_number_MTSD;
    
    %% averages density
    cd(dens_dir);
    Density_file = ['Summary_kmeans_', num2str(Number),'.csv'];
    summary_dens = csvread(Density_file,1,0);
    outlier_area = isoutlier(summary_dens(:,10), 'median');
    outlier_ecc = outlier_area + isoutlier(summary_dens(:,11), 'median');
    outlier_number_dens = length(outlier_ecc(outlier_ecc ~= 0));
    summary_dens(outlier_ecc ~= 0,:) = [];
    summary_dens(isnan(summary_dens(:,3)) == 1,:) = [];
    Averages_dens(loop,1) = Number;
    % Signal area
    Averages_dens(loop,2) = mean(summary_dens(:,2));
    Averages_dens(loop,3) = sqrt(var(summary_dens(:,2))/length(summary_dens(:,2)));
    % Density
    Averages_dens(loop,4) = mean(summary_dens(:,3));
    Averages_dens(loop,5) = sqrt(var(summary_dens(:,3))/length(summary_dens(:,3)));
    % Bundling
    Averages_dens(loop,6) = mean(summary_dens(:,4));
    Averages_dens(loop,7) = sqrt(var(summary_dens(:,4))/length(summary_dens(:,4)));
    % Uniformity
    Averages_dens(loop,8) = mean(summary_dens(:,5));
    Averages_dens(loop,9) = sqrt(var(summary_dens(:,5))/length(summary_dens(:,5)));
    % Sparseness
    Averages_dens(loop,10) = mean(summary_dens(:,6));
    Averages_dens(loop,11) = sqrt(var(summary_dens(:,6))/length(summary_dens(:,6)));
    % Skewness
    Averages_dens(loop,12) = mean(summary_dens(:, 7));
    Averages_dens(loop,13) = sqrt(var(summary_dens(:, 7))/length(summary_dens(:, 7)));
    % Kurtosis
    Averages_dens(loop,14) = mean(summary_dens(:, 8));
    Averages_dens(loop,15) = sqrt(var(summary_dens(:, 8))/length(summary_dens(:, 8)));
    
    
    % gaps
    Averages_dens(loop,16) = mean(summary_dens(:,9));
    Averages_dens(loop,17) = sqrt(var(summary_dens(:,9))/length(summary_dens(:,9)));
    % cell number and outliers
    Averages_dens(loop,18) = length(summary_dens(:,1));
    Averages_dens(loop,19) = outlier_number_dens;
    % area
    Averages_dens(loop,20) = mean(summary_dens(:,10));
    Averages_dens(loop,21) = sqrt(var(summary_dens(:,10))/length(summary_dens(:,10)));
    % eccentricity
    Averages_dens(loop,22) = mean(summary_dens(:,11));
    Averages_dens(loop,23) = sqrt(var(summary_dens(:,11))/length(summary_dens(:,11)));
end

cd(sum_dir);
headers_MTSD = {'Embryo', 'Area','sem',  'Eccentricity','sem', 'Direction_cell','sem',...
    'SD', 'sem', 'Direction_cytoskeleton','sem', 'Aspect ratio','sem', 'Alignment','sem',...
    'Cell number','Outliers'};
csvwrite_with_headers('Summary_MTSD.csv',Averages_MTSD,headers_MTSD);

headers_dens = {'Embryo', 'Signal area', 'sem','Density','sem','Bundling','sem',...
    'Uniformity','sem', 'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'Gaps', 'sem', 'Cell number','Outliers', 'Area', 'sem', 'Eccentricity','sem'};
csvwrite_with_headers('Summary_density.csv',Averages_dens,headers_dens);

%% Writing down pulled data MTSD
Averages_MTSDall = zeros(1,2);
data_k = zeros(1,2);
Averages_MTSDbinned = zeros(1,5);
cd(MTSD_dir);
for loop=1:length(files)
    MTSD_file = ['Summary_MTSD',num2str(loop),'.csv'];
    Data = csvread(MTSD_file,1,0);
    temp = [Data(:,3), Data(:,5)];
    Data2 = Averages_MTSDall;
    Averages_MTSDall  = [Data2; temp];
end
Averages_MTSDall(1,:) = [];
Averages_MTSDall = sortrows(Averages_MTSDall,1);
cd(sum_dir);
csvwrite('MTSDall.csv', Averages_MTSDall);

for k=0.92:0.005:0.99
    data_k =  Averages_MTSDall( Averages_MTSDall(:,1) >= k-0.0025 & Averages_MTSDall(:,1) < k+0.0025,:);
    N = int8((k-0.915)/0.005);
    Averages_MTSDbinned(N,1) = 0;
    Averages_MTSDbinned(N,2) = 0;
    Averages_MTSDbinned(N,3) = 0;
    Averages_MTSDbinned(N,4) = 0;
    if size(data_k,1)>1
        Averages_MTSDbinned(N,1) = mean(data_k(:,1));
        Averages_MTSDbinned(N,2) = sqrt(var(data_k(:,1))/length(data_k(:,1)));
        Averages_MTSDbinned(N,3) = mean(data_k(:,2));
        Averages_MTSDbinned(N,4) = sqrt(var(data_k(:,2))/length(data_k(:,2)));
    elseif size(data_k,1)==1
        Averages_MTSDbinned(N,1) = data_k(:,1);
        Averages_MTSDbinned(N,2) = 0;
        Averages_MTSDbinned(N,3) = data_k(:,2);
        Averages_MTSDbinned(N,4) = 0;
    end
    Averages_MTSDbinned(N,5) = size(data_k,1);
end

headers = {'Eccentricity','sem', 'MTSD', 'sem','cells'};
csvwrite_with_headers('MTSD_binned.csv',Averages_MTSDbinned,headers);


%% writing tissue level MTSD
cd(dist_dir);
MTSDpulled(:,1) = csvread('1_distribution.csv',0,0, [0 0 44 0]);
Eccpulled = zeros(1,1);
counter_cells = 0;
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';
for loop=1:numel(files)
    
    %% reading files
    
    clear Name Number Actin_file Image_actin Path  Image_borders
    cd(image_dir2);
    Number = loop;
    Actin_file = [num2str(Number),'.tif'];
    Image = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(Number),'/'];
    cd(Path);
    Image_borders = imread('tracked_bd.png');
    [im_x, im_y] = size(Image);
    
    %% Collect data about cells and boundaries
    borders;
    
    image_original_double_MTSD = double(im2uint16(Image));
    
    for k = 1:numel(b_valid)
        
        %Apply Sobel Filter over a MTs image to test it
        clear H_full V_full H V M D x y mxd_thr mxd_corrected mxd_indexed
        object_double = image_original_double_MTSD .*...
            poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        H_full = conv2(object_double,Gx);
        V_full = conv2(object_double,Gy);
        H = H_full(5:im_x,5:im_y);
        V = V_full(5:im_x,5:im_y);
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
        mxd_corrected = mxd_corrected + 90 - cell_data(k,5);
        mxd_corrected(mxd_corrected>90,1) = mxd_corrected(mxd_corrected>90,1)-180;
        mxd_corrected(mxd_corrected<-90,1) = mxd_corrected(mxd_corrected<-90,1)+180;
        
        % Make histogram
        [N, bins] = histc(mxd_corrected(:,1),binrange);
        
        mxd_indexed(:,2) = mxd_corrected(:,2);
        mxd_indexed(:,1) = bins(:,1);
        
        
        % Make distribution
        m_added = zeros(45,1);
        for i=1:(length(binrange)-1)
            m_added(i,1) = sum(mxd_indexed(mxd_indexed(:,1) == i,2));
        end
        counter_cells = counter_cells + 1;
        MTSDpulled(:,counter_cells+1) = m_added/sum(m_added);
        Eccpulled(counter_cells) = cell_data(k,4);
    end
end
m_added_norm = zeros(45,16);
m_added_norm(:,1) = MTSDpulled(:,1);
Binnednumber = zeros(15,1);
Binecc = zeros(15,3);
for k=0.92:0.005:0.99
    N = int8((k-0.915)/0.005);
    Binecc_temp = zeros(1,1);
    for i = 1:length(Eccpulled)
       if  Eccpulled(i)>k-0.0025 && Eccpulled(i)<k+0.0025
           m_added_norm(:,N+1) = m_added_norm(:,N+1) + MTSDpulled(:,i+1);
           Binnednumber(N,1) = Binnednumber(N,1) + 1;
           Binecc_temp(Binnednumber(N,1),1) = Eccpulled(i);
       end
    end
    Binecc(N,1) = mean(Binecc_temp(:,1));
    Binecc(N,2) = sqrt(var(Binecc_temp(:,1))/Binnednumber(N,1));
    Binecc(N,4) = Binnednumber(N,1);
    m_added_norm(:,N+1) = m_added_norm(:,N+1)/Binnednumber(N,1);
end

vonmises_fit_dist_sum;
Binecc(:,3) = SD';
Binecc(isnan(Binecc(:,3)) == 1,:) = [];

cd(sum_dir);
headers = {'Eccentricity','sem','MTSD', 'N'};
csvwrite_with_headers('MTSDpulled.csv',Binecc,headers);



cd(currdir);
clc
clear variables
close all