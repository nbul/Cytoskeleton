Change = struct([]);
Cange_norm = zeros(1,1);
for k = 1:numel(b_valid);
    object_double_density = image_original_double .*...
        im2double(poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y));
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    p=0;

    for j = 2:(im_y-1);
        for i=2:(im_x-1);
            %Only directions different to zero are added to table
            if (object_double_density(i,j) && object_double_density(i+1,j) &&...
                    object_double_density(i-1,j) && object_double_density(i,j+1) &&...
                    object_double_density(i,j-1)  && object_double_density(i-1,j-1) &&...
                    object_double_density(i-1 , j+1) && object_double_density(i+1 , j-1) &&...
                    object_double_density(i+1, j+1) ~= 0);                
                p = p + 1;
                Change{k}(p) = (abs(object_double_density(i+1,j) - object_double_density(i,j)) + ...
                    abs(object_double_density(i,j+1) - object_double_density(i,j)) + ...
                    abs(object_double_density(i-1,j) - object_double_density(i,j)) + ...
                    abs(object_double_density(i,j-1) - object_double_density(i,j)) + ...
                    abs(object_double_density(i-1,j-1) - object_double_density(i,j)) + ...
                    abs(object_double_density(i+1,j-1) - object_double_density(i,j)) + ...
                    abs(object_double_density(i-1,j+1) - object_double_density(i,j)) + ...
                    abs(object_double_density(i+1,j+1) - object_double_density(i,j)))/8;
            end
        end;
    end;
    Change{k} =  Change{k}/mean(to_analyse_back_o.PixelValues);
    Cange_norm(k) = length(Change{k}(Change{k}>1))/length(Change{k});
end

edge(Image)
imshow(imoverlay(imadjust(Image2),edge(imadjust(Image2), 'Sobel', 0.05)))


%%Old
to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_all = struct([]);
to_analyse_back_o = struct([]);

mts_density = zeros(1, numel(b_valid));
mts_area = zeros(1, numel(b_valid));
Uniformity = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));
mts_bundling = zeros(1, numel(b_valid));

%% Thresholding an image
image_original_double = im2double(Image2); % original into double
image_edges = edge(Image2, 'Canny');
signal_edges = image_original_double .* image_edges;
[parmhat,parmci] = lognfit(signal_edges(signal_edges>0));
threshold = exp(parmhat(1)-parmhat(2)*parmhat(2));
%threshold2 = adaptthresh(Image2);
im_bin_c = imbinarize(image_original_double,threshold*0.7);
%levels_density = (adaptthresh(image_adjusted)); % treshold value
%levels_density = (graythresh(image_adjusted))*0.7;

%% Generate Masks.
 
signal_original = image_original_double .* im_bin_c;
background_original = image_original_double .* (ones(im_x,im_y) - im_bin_c);

for k = 1:numel(b_valid);
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
    to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    to_analyse_all = regionprops(selected_signal, image_original_double,'PixelValues');
    
    %Cell-by-cell Density Analysis Values
    
    sum_pixvalues_o = sum(to_analyse_o.PixelValues);
    sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
    num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
    num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
    
    mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
    
    if max(to_analyse_c.PixelValues)~= 0
        mts_bundling(k) = mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))...
            /mode(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
    else
        mts_bundling(k) = 0;
    end
    
    Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
        mean(to_analyse_all.PixelValues))/length(to_analyse_all.PixelValues));
   
    Spars(k) = calcSparseness(to_analyse_all.PixelValues,1);
end

%min(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))

%% Skewness
to_analyse_all = struct([]);

Uniformity = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));
cov= zeros(1, numel(b_valid));
cov2= zeros(1, numel(b_valid));
skew= zeros(1, numel(b_valid));

%% Thresholding an image
image_original_double = im2double(Image2); % original into double

%% Generate Masks.
 
for k = 1:numel(b_valid);
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_all = regionprops(selected_signal, image_original_double,'PixelValues');
    
    %Cell-by-cell Density Analysis Values
    [parmhat,parmci] = lognfit(to_analyse_all.PixelValues(to_analyse_all.PixelValues>0));
    cov(k) = sqrt(exp(parmhat(2)*parmhat(2))-1);
    cov2(k) = sqrt((exp(parmhat(2)*parmhat(2))-1)*(exp(2*parmhat(1)+parmhat(2)*parmhat(2))))/...
        exp(parmhat(1));
    skew(k) = (exp(parmhat(2)*parmhat(2))+1)*sqrt(exp(parmhat(2)*parmhat(2))-1);
    
    Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
        mean(to_analyse_all.PixelValues))/length(to_analyse_all.PixelValues));
   
    Spars(k) = calcSparseness(to_analyse_all.PixelValues,1);
end

%min(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))