to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_back_o = struct([]);

mts_density = zeros(1, numel(b_valid));
mts_area = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));

for k = 1:numel(b_valid);
    clear selected_signal
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
    to_analyse_c = regionprops(selected_signal, signal_corrected,'PixelValues');
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    
    %Cell-by-cell Density Analysis Values
    
    sum_pixvalues_o = sum(to_analyse_o.PixelValues(:,1));
    sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
    num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
    num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
    
    mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
    
    Spars(k) = calcSparseness(to_analyse_c.PixelValues,1);
end