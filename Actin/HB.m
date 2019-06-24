phiarray = 2:4:178;

bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
m_added_norm = zeros(45,length(0.7:0.01:0.98)+1);
m_added_norm(:,1) = bincenter;

for i = 1:1:29
    ecc = 0.69 + 0.01 * i;
    a = 1; b = a/sqrt(1-ecc^2);
    fHB = 1./ (cosd(phiarray).^2 + b^2 * sind(phiarray).^2);
    fHB_norm = fHB./trapz(phiarray,fHB );
    m_added_norm(:,i+1) = fHB_norm';
    
end
vonmises;
result = zeros(29,2);
result(:,1) = 0.69 + 0.01 * (1:1:29);
result(:,2) = SD';