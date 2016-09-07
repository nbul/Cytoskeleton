%% This script update is meant to analyse density with a normalisation step included.

% Collect data about cells and boundaries

clear image_borders checknumbers borders_bin B L s
borders_bin = im2bw(Image_borders,0);
[B,L] = bwboundaries(borders_bin,'holes');
figure, imshow(label2rgb(L, @jet, [.5 .5 .5]));
hold on;

s = regionprops(L, 'Centroid');
for k = 1:length(B);
    clear boundary
    boundary = B{k};
    clear c
    c = s(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2);
end;

clear im_cells_data
im_cells_data=regionprops(L,'Area', 'Eccentricity','MajorAxisLength','MinorAxisLength' , 'Orientation','PixelIdxList','PixelList','Centroid');


%% Recording cells and their boundaries for cells larger than 200 px
clear b_valid
cell_counter = 5;
cell_filteredcounter = 0;
cell_data = zeros(numel(im_cells_data), 4);
good_cell = zeros(1, numel(im_cells_data));

b_valid = cell(0);
for i=6:numel(im_cells_data);
    if im_cells_data(i).Area > 200;
        cell_counter = cell_counter + 1;
        cell_filteredcounter = cell_filteredcounter + 1;
        cell_data(cell_filteredcounter,1) = cell_counter;
        cell_data(cell_filteredcounter,2) = im_cells_data(i).Area;
        cell_data(cell_filteredcounter,3) = im_cells_data(i).Eccentricity;
        cell_data(cell_filteredcounter,4) = im_cells_data(i).Orientation;
        good_cell(cell_filteredcounter)=i;
        b_valid(cell_filteredcounter) = B(i);
    end
end
cell_data(cell_data(:,2)==0,:) = [];
good_cell(:,good_cell(1,:)==0) = [];
b_valid = b_valid(~cellfun('isempty',b_valid));  

%% Open MTs image. Adjust image to optimal settings.

imbackground = imopen(Image_actin,strel('disk',40)); %background
imcorrected = Image_actin - imbackground; % background subtraction
imcorrected_double = im2double(imcorrected); % background subtracted into double
image_original_double = im2double(Image_actin); % original into double
image_adjusted = imadjust(imcorrected); %adjusting intensity
image_adjusted_double = im2double(image_adjusted); % with adjusted intensity into double

[im_x, im_y] = size(image_original_double);

%% Processed Image for Density Analysis

levels_density = (graythresh(image_adjusted))*1; % treshold value
im_bin_c = im2bw(image_adjusted,levels_density); % thresholding adjusted image

%% Drowing overlay of selected cell borders on image of cytoskeleton
image1=figure;
imshow(image_adjusted), title('Adjusted MTs Image');
hold on;
for k = 1:length(b_valid);
    clear boundary_valid
    boundary_valid = b_valid{k};
    
    clear c
    c = im_cells_data(good_cell(k)).Centroid;
    
    clear c_labels
    c_labels = text(c(1), c(2), sprintf('%d', k),'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    set(c_labels,'Color',[1 1 0])
    plot(boundary_valid(:,2), boundary_valid(:,1), 'r', 'LineWidth', 2);
end;

%% Generate Cell Masks.

clear select_object1 object_double mts_area to_analyse to_analyse_o to_analyse_c
clear to_analyse_back_o to_analyse_back_c
clear signal_original signal_corrected
clear background_im_bin_c background_original background_corrected
clear mts_density result m_added_norm

inv_one_matrix = ones(im_x,im_y);
background_im_bin_c = inv_one_matrix - im_bin_c;
signal_original = image_original_double .* im_bin_c;
signal_corrected = (im2double(imcorrected)) .* im_bin_c;
background_original = image_original_double .* background_im_bin_c;
background_corrected = (im2double(imcorrected)) .* background_im_bin_c;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';

%%  Assign memory.
select_object1 = struct([]);
to_analyse_o = struct([]);
to_analyse_c = struct([]);
to_analyse_back_o = struct([]);
to_analyse_back_c= struct([]);
result = struct([]);

m_added_norm = zeros(45,numel(b_valid)+1);
mts_density = zeros(1, numel(b_valid));
mts_area = zeros(1, numel(b_valid));
%%  Conduct cell-by-cell analysis.
for k = 1:numel(b_valid);
    
    select_object1{k} = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    object_double = im2double(select_object1{k});
    
    % Density data from signal and background
    
    to_analyse_o{k} = regionprops(select_object1{k},signal_original,'PixelList','PixelValues','Area');
    to_analyse_c{k} = regionprops(select_object1{k},signal_corrected,'PixelList','PixelValues','Area');
    to_analyse_back_o{k} = regionprops(select_object1{k},background_original,'PixelList','PixelValues','Area');
    to_analyse_back_c{k} = regionprops(select_object1{k},background_corrected,'PixelList','PixelValues','Area');
    
    %Cell-by-cell Density Analysis Values
    
    sum_pixvalues_o = sum(to_analyse_o{k}.PixelValues(:,1));
    sum_pixvalues_back_o = sum(to_analyse_back_o{k}.PixelValues(:,1));
    num_pixvalues_c = 0;
    num_pixvalues_back_c = 0;
    
    for density_values = 1 : numel(to_analyse_o{k}.PixelValues(:,1));
        if to_analyse_c{k}.PixelValues(density_values,1) ~= 0;
            num_pixvalues_c = num_pixvalues_c + 1;
        end
        if to_analyse_back_c{k}.PixelValues(density_values,1) ~= 0;
            num_pixvalues_back_c = num_pixvalues_back_c + 1;
        end
    end
    mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
    
    %Apply Sobel Filter over a MTs image to test it
    clear H_full V_full H V M D_radians pp qq x y
    result{k} = image_original_double .* object_double;
    H_full = conv2(result{k},Gx);
    V_full = conv2(result{k},Gy);
    H = H_full(5:im_x,5:im_y);
    V = V_full(5:im_x,5:im_y);
    M = sqrt(H.^2 + V.^2);
    
    clear D_radians D 
    D_radians = atan2(V, H);
    D = -(180/pi) * D_radians;
    
    [x, y] = size(M);
    p = 1;
    mxd = zeros(im_cells_data(good_cell(k)).Area,2);
    for j = 2:(y-1);
        for i=2:(x-1);
            if ((M(i,j)) & (M(i+1,j)) & (M(i-1,j)) & (M(i,j+1))  & (M(i,j-1) ~=0) & (M(i-1,j-1)) & (M(i-1 , j+1)) & (M(i+1 , j-1)) & (M(i+1, j+1))) ~= 0 ; %Only directions different to zero are added to table
                mxd(p,2) = M(i,j); %Second column with magnitudes
                mxd(p,1) = D(i,j); %First column with angles
                p = p + 1;
            end;
        end;
    end;
    mxd(mxd(:,2)==0,:) = [];
    
    mxd_corrected = zeros(length(mxd),2);

    [r, ~] = size(mxd);
    max_mxd  = max(mxd(:,2)); %maximum magnitude
    mxd_thr = mxd;
    mxd_thr(:,2) = mxd_thr(:,2)/max_mxd; %normalised to max magnitude
    
    % Remove all pixels with magnitude less than 22% of maximum
    p = 1;
    for i = 1 : r;
        if mxd_thr(i,2) >= 0.22;
            mxd_corrected(p,1) = mxd_thr(i,1);
            mxd_corrected(p,2) = mxd_thr(i,2);
            p = p + 1;
        end
    end
    mxd_corrected(mxd_corrected(:,2)==0,:) = [];

    [r, ~] = size(mxd_corrected);
    
    for i = 1 : r;
        if mxd_corrected(i,1) < 0;
            mxd_corrected(i,1) = mxd_corrected(i,1) + 180;
        end;
    end;
    
    % Sorting by angle and shifting to -90 to 90
    mxd_corrected = sortrows(mxd_corrected,1);
    mxd_corrected = [mxd_corrected(:,1) - 90, mxd_corrected(:,2)]; %To place the values between [-90 90]
    
    for xx = 1 : numel(mxd_corrected(:,1));
        if mxd_corrected(xx,1) >= 90;
            mxd_corrected(xx,1) = 89.9;
        end;
    end;
    
    mxd_corrected = sortrows(mxd_corrected,1);
    
    % Make histogram
    
 
    [N, bins] = histc(mxd_corrected(:,1),binrange);
   
    clear mxd_indexed
    mxd_indexed(:,2) = mxd_corrected(:,2);
    mxd_indexed(:,1) = bins(:,1);

    
    % Make distribution
    m_added = zeros(45,2);
    for i=1:(length(binrange)-1)
        m_added(i,1) = bincenter(i);
        for ii = 1:r
            if mxd_indexed(ii,1) == i
                m_added(i,2) = m_added(i,2) + mxd_indexed(ii,2);
            end
        end
    end
    
    total = sum(m_added);
    m_added_norm(:,1) = m_added(:,1);
    m_added_norm(:,k+1) = m_added(:,2)/total(2);
end;


























