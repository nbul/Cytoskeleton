%% Open MTs image. Adjust image to optimal settings.

imcorrected = Image2 - imopen(Image2,strel('disk',40)); % background subtraction
image_original_double = im2double(Image2); % original into double
image_original_double_MTSD = im2double(Image);

%% Processed Image for Density Analysis
image_adjusted = imadjust(imcorrected); %adjusting intensity
levels_density = (graythresh(image_adjusted))*1; % treshold value
im_bin_c = im2bw(image_adjusted,levels_density); % thresholding adjusted image

%% Generate Cell Masks.
signal_original = image_original_double .* im_bin_c;
signal_corrected = (im2double(imcorrected)) .* im_bin_c;
background_original = image_original_double .* (ones(im_x,im_y) - im_bin_c);

%% Drowing overlay of selected cell borders on image of cytoskeleton
image1=figure;
imshow(image_adjusted), title('Adjusted MTs Image');
hold on;
for k = 1:length(b_valid);
    clear boundary_valid
    boundary = b_valid{k};
    c = im_cells_data(cell_data(k,2)).Centroid;
    c_labels = text(c(1), c(2), sprintf('%d', k),'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Fontsize', 10);
    set(c_labels,'Color',[1 1 0])
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
end;

%% Saving the image
cd(im_dir);
image_filename = [num2str(Number),'_analysed_image.tif'];
print(image1, '-dtiff', '-r150', image_filename);
close all



