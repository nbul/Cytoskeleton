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

%% Processed Image for Density Analysis
if method == 1
    image_original_double = double(im2uint16(Image2)); % original into double
else
    image_original_double = im2double(Image2);
end

image_edges = edge(Image2, 'Canny');
signal_edges = image_original_double .* image_edges;

if method == 1
    [A1, A2] = histcounts(signal_edges(signal_edges>0));
    bins = (min(A2)+((A2(2)-A2(1))/2)):((A2(2)-A2(1))):(max(A2)-((A2(2)-A2(1))/2));
    [XOut,YOut] = prepareCurveData(bins,A1);
    fo = fitoptions('gauss1', 'Lower', [0 min(A2) 0], 'Upper', [Inf max(A2) Inf]);
    [threshold, gof_edges] = fit(XOut, YOut, 'gauss1', fo);
    im_bin_c = imbinarize(image_original_double,threshold.b1*0.7);
else
    im_bin_c = imbinarize(imadjust(image_original_double),graythresh(imadjust(image_original_double))*0.7);
end

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
end

