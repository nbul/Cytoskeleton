%% Preallocate memory

to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_all = struct([]);
to_analyse_back_o = struct([]);

mts_area = zeros(1, numel(b_valid));
mts_signal = zeros(1, numel(b_valid));

%% Processed Image for Density Analysis
image_original_double = im2double(Image);
threshold = graythresh(imadjust(image_original_double))*0.7;
im_bin_c = imbinarize(imadjust(image_original_double),graythresh(imadjust(image_original_double))*0.7);

%% Generate Cell Masks.
signal_original = image_original_double .* im_bin_c;
background_original = image_original_double .* (ones(im_x,im_y) - im_bin_c);

for k = 1:numel(b_valid)
    clear selected_signal
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
    to_analyse_o = regionprops(selected_signal, Image,'PixelValues');
    
    % Relative to background and signal area 
    num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
    num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
    
    
    % Signal Area
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);  
    mts_signal(k) = mean(to_analyse_o.PixelValues(to_analyse_o.PixelValues(:, 1) ~= 0,1))/...
        mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
end

