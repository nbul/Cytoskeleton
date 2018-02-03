%% Preallocate memory


to_analyse_all = struct([]);

Uniformity = zeros(1, numel(b_valid));
Spars = zeros(1, numel(b_valid));
kurt = zeros(1, numel(b_valid));
skew = zeros(1, numel(b_valid));
Sdr = zeros(1, numel(b_valid));
Sdq = zeros(1, numel(b_valid));
image_original_double = im2double(Image2);
[px, py] = gradient(image_original_double);

for k = 1:numel(b_valid)
    clear selected_signal
    % Density data from signal and background
    selected_signal = poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
    to_analyse_all = regionprops(selected_signal, image_original_double,'PixelValues');
    signal_original = image_original_double .* selected_signal;
    px2 = px .* selected_signal;
    py2 = py .* selected_signal;
    
    % Uniformity
    Uniformity(k) = 100 * (1 - sqrt(var(to_analyse_all.PixelValues))/mean(to_analyse_all.PixelValues));
    %Sparseness
    Spars(k) = calcSparseness(to_analyse_all.PixelValues/mean(to_analyse_all.PixelValues),1);
    
    % Kurtosis and skewness
    kurt(k) = kurtosis(to_analyse_all.PixelValues);
    skew(k) = skewness(to_analyse_all.PixelValues);
    
    %Sdr
    Sdr(k) = (sqrt(1+sum(px(:).*px(:)+(py(:).*py(:))))-1)/length(to_analyse_all.PixelValues);
    %Sdq - the root mean square gradient
    Sdq(k) = sqrt(sum(px(:).*px(:)+(py(:).*py(:)))/length(to_analyse_all.PixelValues));
    
    
end

