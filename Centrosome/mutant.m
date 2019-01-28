cd(b_dir_mutant);
files_mutant = dir('*.tif');
cd(currdir);
Averages = zeros(length(files_mutant), 16);
Pulled = zeros(1,8);
for loop=1:length(files_mutant)
    
    cd(b_dir_mutant);
    Path = [b_dir_mutant, '/', num2str(loop),'/'];
    cd(Path);
    Image_borders = imread('handCorrection.tif');
    [im_x, im_y] = size(Image);
    %% Border data
    cell_data = zeros(1, 8);
    borders;
    
    %% Centrosome image
    cd(cent_dir_mutant);
    Cs_file = [num2str(loop),'.tif'];
    Image = imread(Cs_file);
    Image = imgaussfilt(Image,2);
    image_original_double = im2double(Image);
    im_bin_c = imbinarize(imadjust(image_original_double),thr);
    im_bin_c = bwareaopen(im_bin_c,10);
    im_bin_c = imfill(im_bin_c,'holes');
    imageprocessing;
    
    %% Saving the image
    cd(im_dir_mutant);
    image_filename = [num2str(loop),'_analysed_image.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    
    for k = 1:numel(b_valid)
        object_double = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        Cs_count_image = object_double .* im_bin_c;
        CC = bwconncomp(Cs_count_image);
        S = regionprops(CC, 'Centroid');
        S2 = regionprops(CC,Image, 'PixelValues');
        if numel(S)>2
            Cs_count_image = bwareafilt(imbinarize(Cs_count_image),2,'largest');
            CC = bwconncomp(Cs_count_image);
            S = regionprops(CC, 'Centroid','Area');
            S2 = regionprops(CC,Image, 'PixelValues');
        end
        cell_data(k,6) = numel(S);
        if numel(S) == 2
            cell_data(k,8) = sqrt((S(1).Centroid(1,1)-S(2).Centroid(1,1))*(S(1).Centroid(1,1)-S(2).Centroid(1,1)) +...
                (S(1).Centroid(1,2)-S(2).Centroid(1,2))*(S(1).Centroid(1,2)-S(2).Centroid(1,2)));
        end
        for m = 1:numel(S)
            cell_data(k,7) = cell_data(k,7) + sum(S2(m).PixelValues);
        end
    end
    
    %% Cell by cell summary
    cd(sum_dir_mutant_cs);
    summary_filename = ['Summary_',num2str(loop),'.csv'];
    csvwrite_with_headers(summary_filename,cell_data,headers2);
    
    %% Average values
    Averages(loop,1) = loop;
    Averages(loop,2:2:size(cell_data(:,2:end),2)*2) = mean(cell_data(:,2:end),1);
    Averages(loop,3:2:(size(cell_data(:,2:end),2)*2+1))...
        = sqrt(var(cell_data(:,2:end),1)/size(cell_data(:,2:end),1));
    Averages(loop,size(cell_data(:,2:end),2)*2+2) = size(cell_data(:,2:end),1);
    
    %% Pulled data
    Pulled = [Pulled;cell_data];
end

Pulled(Pulled(:,1) == 0,:) = [];
cd(sum_dir_mutant);
summary_filename = 'Summary_average.csv';
csvwrite_with_headers(summary_filename,Averages,headers1);
summary_filename = 'Summary_pulled.csv';
csvwrite_with_headers(summary_filename,Pulled,headers2);
cd(currdir);