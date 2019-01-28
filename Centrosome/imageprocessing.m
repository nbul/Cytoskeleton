%% Drowing overlay of selected cell borders on image of cytoskeleton
image1=figure;
imshow(im_bin_c);
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





