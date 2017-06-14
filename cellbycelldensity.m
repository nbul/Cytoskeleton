to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_back_o = struct([]);
to_analyse_back_c= struct([]);

mts_density = zeros(1, numel(b_valid));
mts_area = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));

for k = 1:numel(b_valid);
    
    % Density data from signal and background
    
    to_analyse_o = regionprops(poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y),...
        signal_original,'PixelValues');
    to_analyse_c = regionprops(poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y),...
        signal_corrected,'PixelValues');
    to_analyse_back_o = regionprops(poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y),...
        background_original,'PixelValues');
    to_analyse_back_c = regionprops(poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y),...
        background_corrected,'PixelValues');
    
    %Cell-by-cell Density Analysis Values
    
    sum_pixvalues_o = sum(to_analyse_o.PixelValues(:,1));
    sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
    num_pixvalues_c = 0;
    num_pixvalues_back_c = 0;
    
    for density_values = 1 : numel(to_analyse_o.PixelValues(:,1));
        if to_analyse_c.PixelValues(density_values,1) ~= 0;
            num_pixvalues_c = num_pixvalues_c + 1;
        end
        if to_analyse_back_c.PixelValues(density_values,1) ~= 0;
            num_pixvalues_back_c = num_pixvalues_back_c + 1;
        end
    end
    mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
    
    to_analyse_c.PixelValues(to_analyse_c.PixelValues==0)=[];
    Spars(k) = calcSparseness(to_analyse_c.PixelValues,1);
end