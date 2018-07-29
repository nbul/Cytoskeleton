clear object_double  
image_original_double_MTSD = double(im2uint16(Image));
%%  Assign memory.

m_added_norm = zeros(45,numel(b_valid)+1);
m_added_norm(:,1) = bincenter;

%%  Conduct cell-by-cell analysis.
for k = 1:numel(b_valid)
    
    %Apply Sobel Filter over a MTs image to test it
    clear H_full V_full H V M D x y mxd_thr mxd_corrected mxd_indexed
    object_double = image_original_double_MTSD .*...
        poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    H_full = conv2(object_double,Gx);
    V_full = conv2(object_double,Gy);
    H = H_full(5:im_x,5:im_y);
    V = V_full(5:im_x,5:im_y);
    M = sqrt(H.^2 + V.^2);
    D = -(180/pi) * atan2(V, H);
    
    [x, y] = size(M);
    p = 1;
    mxd = zeros(1,2);
    
    for j = 2:(y-1)
        for i=2:(x-1)
            %Only directions different to zero are added to table
            if ((M(i,j)) & (M(i+1,j)) & (M(i-1,j)) & (M(i,j+1))  & (M(i,j-1) ~=0) &...
                    (M(i-1,j-1)) & (M(i-1 , j+1)) & (M(i+1 , j-1)) & (M(i+1, j+1))) ~= 0  
                mxd(p,2) = M(i,j); %Second column with magnitudes
                mxd(p,1) = D(i,j); %First column with angles
                p = p + 1;
            end
        end
    end
    max_mxd  = max(mxd(:,2)); %maximum magnitude
    mxd_thr = mxd./repmat([1,max_mxd], length(mxd), 1); %normalised to max magnitude
    
    % Remove all pixels with magnitude less than 22% of maximum
    
    mxd_corrected = mxd_thr(mxd_thr(:,2) >= 0.22,:);
    mxd_corrected(mxd_corrected(:,1) < 0,1) = mxd_corrected(mxd_corrected(:,1) < 0,1) + 180;
    mxd_corrected(:,1) = mxd_corrected(:,1) - 90;
    mxd_corrected(mxd_corrected(:,1) >= 90,1) = 89.9;
    mxd_corrected = sortrows(mxd_corrected,1);
    
    % Make histogram
    [N, bins] = histc(mxd_corrected(:,1),binrange);
    
    mxd_indexed(:,2) = mxd_corrected(:,2);
    mxd_indexed(:,1) = bins(:,1);
    
    
    % Make distribution
    m_added = zeros(45,1);
    for i=1:(length(binrange)-1)
        m_added(i,1) = sum(mxd_indexed(mxd_indexed(:,1) == i,2));
    end
   
    m_added_norm(:,k+1) = m_added/sum(m_added);
end

%% writing down distributions in each image
cd(dist_dir);
gradient_filename = [num2str(Number),'_distribution.csv'];
csvwrite(gradient_filename, m_added_norm);

























