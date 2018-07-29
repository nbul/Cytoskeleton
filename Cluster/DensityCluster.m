%% Processed Image for Density Analysis
image_original_double = im2double(Image2);

ImageX = Image2(:);

parfor cl = 2:20
    clust2(:,cl-1) = kmeans(im2double(ImageX), cl, 'replicate',5);
end
eva = evalclusters(im2double(ImageX),clust2,'DaviesBouldin');
km = kmeans(im2double(ImageX),eva.OptimalK,'replicate',5);
km2 = reshape(km, im_y,im_x);

thr = zeros(eva.OptimalK,1);
for clust = 1:eva.OptimalK
    thr(clust) = mean(Image2(km2==clust));
end
[Num1, Idx1] = min(thr);
threshold = max(Image2(km2==Idx1));
im_bin_c = imbinarize(im2double(Image2),double(threshold)/255/255);
im_bin_b = imcomplement(im_bin_c);

data_dens = [eva.OptimalK, threshold, max(Image2(:))];
headers4 = {'number clusters', 'threshold', 'max intensity'};
cd(dens_dir);
clustering_filename = ['Clustering_',num2str(Number),'.csv'];
csvwrite_with_headers(clustering_filename,data_dens,headers4);
cd(currdir);


