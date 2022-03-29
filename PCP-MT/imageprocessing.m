%% Drowing overlay of selected cell borders on image of cytoskeleton
image1=figure;
imshow(imadjust(Image)), title('Adjusted MTs Image');
hold on;
for k = 1:length(b_valid)
    clear boundary_valid
    boundary = b_valid{k};
    c = im_cells_data(cell_data(k,2)).Centroid;
    c_labels = text(c(1), c(2), sprintf('%d', k),'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Fontsize', 10);
    set(c_labels,'Color',[1 1 0])
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
end

%% Saving the image
cd(im_dir);
image_filename = [num2str(loop),'_analysed_image.tif'];
print(image1, '-dtiff', '-r150', image_filename);
close all


image2=figure;
imshow(imadjust(Image)), title('Adjusted MTs Image');
hold on;
MA = mean([im_cells_data.MajorAxisLength]);
for k = 1:length(b_valid)
    clear boundary_valid
    boundary = b_valid{k};
    c = im_cells_data(cell_data(k,2)).Centroid;
    line([c(1) + cosd(mu(k))*MA/4, c(1) - cosd(mu(k))*MA/4],...
        [c(2) - sind(mu(k))*MA/4, c(2) + sind(mu(k))*MA/4],...
        'Color', 'g', 'LineWidth', 3);
    plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2);
end

%% Saving the image
cd(im_dir2);
image_filename = [num2str(loop),'_average_MT.tif'];
print(image2, '-dtiff', '-r150', image_filename);
close all

image3=figure;
imshow(Image * 0), title('Adjusted MTs Image');
set(gcf, 'InvertHardCopy', 'off');
hold on;
for k = 1:length(b_valid)
    clear boundary_valid
    boundary = b_valid{k};
    mu2(k) = im_cells_data(cell_data(k,2)).Orientation;
    c = im_cells_data(cell_data(k,2)).Centroid;
    line([c(1) + cosd(mu2(k))*MA/4, c(1) - cosd(mu2(k))*MA/4],...
        [c(2) - sind(mu2(k))*MA/4, c(2) + sind(mu2(k))*MA/4],...
        'Color', 'w', 'LineWidth', 3);
    plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2);
end

%% Saving the image
cd(im_dir3);
image_filename = [num2str(loop),'_average_CellOr.tif'];
print(image3, '-dtiff', '-r150', image_filename);
close all

image4=figure;
imshow(Image * 0), title('Adjusted MTs Image');
set(gcf, 'InvertHardCopy', 'off');
hold on;
for k = 1:length(b_valid)
    clear boundary_valid
    boundary = b_valid{k};
    mu2(k) = im_cells_data(cell_data(k,2)).Orientation;
    c = im_cells_data(cell_data(k,2)).Centroid;
    line([c(1) + cosd(mu(k))*MA/4, c(1) - cosd(mu(k))*MA/4],...
        [c(2) - sind(mu(k))*MA/4, c(2) + sind(mu(k))*MA/4],...
        'Color', 'g', 'LineWidth', 3);
    line([c(1) + cosd(mu2(k))*MA/4, c(1) - cosd(mu2(k))*MA/4],...
        [c(2) - sind(mu2(k))*MA/4, c(2) + sind(mu2(k))*MA/4],...
        'Color', 'w', 'LineWidth', 3);
    plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2);
    
end

%% Saving the image
cd(im_dir3);
image_filename = [num2str(loop),'_average_CellMTOr.tif'];
print(image4, '-dtiff', '-r150', image_filename);
close all

if choice == 0
    image5=figure;
    imshow(Image * 0), title('Adjusted MTs Image');
    set(gcf, 'InvertHardCopy', 'off');
    hold on;
    for k = 1:length(b_valid)
        clear boundary_valid
        boundary = b_valid{k};
        c = im_cells_data(cell_data(k,2)).Centroid;
        line([c(1) + cosd(PCPangle(k))*MA/4, c(1) - cosd(PCPangle(k))*MA/4],...
            [c(2) - sind(PCPangle(k))*MA/4, c(2) + sind(PCPangle(k))*MA/4],...
            'Color', '#00FFFF', 'LineWidth', 3);
        plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2);
    end
    
    %% Saving the image
    cd(im_dir4);
    image_filename = [num2str(loop),'_average_PCPOr.tif'];
    print(image5, '-dtiff', '-r150', image_filename);
    close all
    
    image6=figure;
    imshow(Image * 0), title('Adjusted MTs Image');
    set(gcf, 'InvertHardCopy', 'off');
    hold on;
    for k = 1:length(b_valid)
        clear boundary_valid
        boundary = b_valid{k};
        mu2(k) = im_cells_data(cell_data(k,2)).Orientation;
        c = im_cells_data(cell_data(k,2)).Centroid;
        line([c(1) + cosd(PCPangle(k))*MA/4, c(1) - cosd(PCPangle(k))*MA/4],...
            [c(2) - sind(PCPangle(k))*MA/4, c(2) + sind(PCPangle(k))*MA/4],...
            'Color', '#00FFFF', 'LineWidth', 3);
        line([c(1) + cosd(mu2(k))*MA/4, c(1) - cosd(mu2(k))*MA/4],...
            [c(2) - sind(mu2(k))*MA/4, c(2) + sind(mu2(k))*MA/4],...
            'Color', 'w', 'LineWidth', 3);
        plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2);
        
    end
    
    %% Saving the image
    cd(im_dir4);
    image_filename = [num2str(loop),'_average_CellPCPOr.tif'];
    print(image6, '-dtiff', '-r150', image_filename);
    close all
end