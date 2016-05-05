%hists is a 256x256 2D-histogram of grayscale value and neighborhood average grayscale value pair.
%total is the number of pairs in the given image.
%threshold is the threshold obtained.

function threshold = otsu2D(hists, total)
maximum = 0.0;
threshold = 0;
w_0=0;mu_ti=0;mu_tj=0;mu_i=0;mu_j=0;mu_is=0;mu_js=0;
for i =1:256
    for j=1:256
        P(i,j)=hists(i,j)/total;
    end
end
for i=1:256
    for j=1:256
        mu_ti=mu_ti+i*P(i,j);
        mu_tj=mu_tj+j*P(i,j);
    end
end
for ii=1:256
    for jj=1:256
        w_0=w_0+P(ii,jj);
        mu_is=mu_is+ii*P(ii,jj);
        mu_i(ii,jj)=mu_is;
        mu_js=mu_js+jj*P(ii,jj);
        mu_j(ii,jj)=mu_js;
        tr = ((mu_i(ii,jj)-(w_0*mu_ti))^2 + (mu_j(ii,jj)-(w_0*mu_tj))^2)/(w_0*(1-w_0));

        if ( tr >= maximum )
            threshold = ii;
            maximum = tr;
        end
    end
end