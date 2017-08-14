function threshold = TwoDOtsu(hists, total)
threshold = 0;
tr = zeros(256,256);
w0 = 0;
w = zeros(256,256);
m0i = 0;
mi = zeros(256,256);
m0j = 0;
mj = zeros(256,256);
mTi = 0;
mTj = 0;
for ii = 1:256
    for jj = 1:256
        mTi = mTi + ii * hists(ii,jj)/total;
        mTj = mTj + jj * hists(ii,jj)/total;
    end
end

for ii = 1:256
    for jj = 1:256
        w0 = w0 + hists(ii,jj)/total;
        w(ii,jj) = w0;
        m0i = m0i + ii * hists(ii,jj)/total;
        mi(ii,jj) = m0i;
        m0j = m0j + jj * hists(ii,jj)/total;
        mj(ii,jj) = m0j;
        tr(ii,jj) = ((mTi * w0 - m0i)^2 + (mTj * w0 - m0j)^2)/(w0*(1-w0));
    end
end
[M,I] = max(tr(:));
if M>0
    [threshold, I_col] = ind2sub(size(tr),I);
end
end