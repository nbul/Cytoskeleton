%% Preallocate memory

to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_all = struct([]);
to_analyse_back_o = struct([]);

mts_density = zeros(1, numel(b_valid));
mts_area = zeros(1, numel(b_valid));
Uniformity = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));
mts_bundling = zeros(1, numel(b_valid));
kurt = zeros(1, numel(b_valid));
skew = zeros(1, numel(b_valid));
space = zeros(1, numel(b_valid));

%% Processed Image for Density Analysis
image_original_double = im2double(Image2);

ImageX = Image2(:);

parfor cl = 2:20
    clust2(:,cl-1) = kmeans(im2double(ImageX), cl, 'replicate',5);
end
eva = evalclusters(im2double(ImageX),clust2,'DaviesBouldin');
km = kmeans(im2double(ImageX),eva.OptimalK,'replicate',5);
km2 = reshape(km, im_y,im_x);

thr = zeros(eva.OptimalK,1);
for clust = 1:eva.OptimalK
    thr(clust) = mean(Image2(km2==clust));
end
[Num1, Idx1] = min(thr);
threshold = max(Image2(km2==Idx1));
im_bin_c = imbinarize(im2double(Image2),double(threshold)/255/255);
im_bin_b = imcomplement(im_bin_c);

data_dens = [eva.OptimalK, threshold, max(Image2(:))];
headers4 = {'number clusters', 'threshold', 'max intensity'};
cd(dens_dir);
clustering_filename = ['Clustering_',num2str(Number),'.csv'];
csvwrite_with_headers(clustering_filename,data_dens,headers4);
cd(currdir);

%% Generate Cell Masks.
signal_original = image_original_double .* im_bin_c;
background_original = image_original_double .* (ones(im_x,im_y) - im_bin_c);

for k = 1:numel(b_valid)
    clear selected_signal
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
    to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    to_analyse_all = regionprops(selected_signal, image_original_double,'PixelValues');
    
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
            min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))) / ...
            min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
%         mts_bundling(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
%             (sum_pixvalues_back_o / num_pixvalues_back_c));
        %         mts_bundling(k) = mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))...
        %             /min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
    else
        mts_bundling(k) = 0;
    end
    
    % Uniformity
    Uniformity(k) = 100 * (1 - sqrt(var(to_analyse_all.PixelValues))/mean(to_analyse_all.PixelValues));
    %Sparseness
    Spars(k) = calcSparseness(to_analyse_o.PixelValues/mean(to_analyse_o.PixelValues(to_analyse_o.PixelValues>0)),1);
    
    % Kurtosis and skewness
    if max(to_analyse_c.PixelValues)~= 0
        signal = to_analyse_all.PixelValues - min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
        kurt(k) = kurtosis(signal);
        skew(k) = skewness(signal);
    else
        signal = 0;
        kurt(k) = 0;
        skew(k) = 0;
    end
    
    % Spaces
    ccbg = bwconncomp(selected_signal.*im_bin_b);
    Space1 = regionprops(ccbg, 'Area');
    if isempty(Space1)==0
        space(k) = mean(cat(1, Space1.Area));
    else
        space(k) =0;
    end
end

