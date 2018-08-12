%% Preallocate memory

to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_all = struct([]);
to_analyse_back_o = struct([]);

mts_density = zeros(1, numel(b_valid));
mts_area = zeros(1, numel(b_valid));
Uniformity = zeros(1, numel(b_valid));
UNAAD = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));
mts_bundling = zeros(1, numel(b_valid));
kurt = zeros(1, numel(b_valid));
skew = zeros(1, numel(b_valid));
space = zeros(1, numel(b_valid));
Ent1 = zeros(1, numel(b_valid));
Ent2 = zeros(1, numel(b_valid));
Sdq = zeros(1, numel(b_valid));
Sdr = zeros(1, numel(b_valid));
intensity = zeros(1, numel(b_valid));
intensity_mts = zeros(1, numel(b_valid));

%% Processed Image for Density Analysis
im_bin_c = imbinarize(im2double(Image2),double(threshold)/255/255);
im_bin_b = imcomplement(im_bin_c);
cd(currdir);

%% Generate Cell Masks.
signal_original = im2double(Image2) .* im_bin_c;
background_original = im2double(Image2) .* (ones(im_x,im_y) - im_bin_c);

for k = 1:numel(b_valid)
    clear selected_signal
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
    to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    to_analyse_all = regionprops(selected_signal, im2double(Image2),'PixelValues');
    
    % Relative to background and signal area 
    sum_pixvalues_o = sum(to_analyse_o.PixelValues(:,1));
    sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
    num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
    num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
    
    intensity(k) = mean(to_analyse_o.PixelValues)*255*255;
    intensity_mts(k) = mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))*255*255;
    
    %density
    if num_pixvalues_back_c ~= 0
        mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
            (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    else
        mts_density(k) = 255*255*(sum_pixvalues_o / num_pixvalues_c - threshold/255/255)/ threshold;
    end
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
    Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
        (to_analyse_all.PixelValues + mean(to_analyse_all.PixelValues)))/length(to_analyse_all.PixelValues));
    UNAAD(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))/...
        mean(to_analyse_all.PixelValues))/length(to_analyse_all.PixelValues));
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
    
    if max(to_analyse_c.PixelValues)~= 0
        Ent1(k) = entropy(im2double(Image2(selected_signal~= 0)));
        Ent2(k) = entropy(im2double(Image2((selected_signal.*im_bin_c)~= 0)));
    else
        Ent1(k) = 0;
        Ent2(k) = 0;
    end
    
    %Sdq - the root mean square gradient
    [px, py] = gradient(im2double(Image2));
    Sdq(k) = sqrt(sum(px(:).*px(:)+(py(:).*py(:)))/length(to_analyse_all.PixelValues));
    %Sdr - the developed interfacial area ratio    
    Sdr(k) = (sqrt(1+sum(px(:).*px(:)+(py(:).*py(:))))-1)/length(to_analyse_all.PixelValues);
end

