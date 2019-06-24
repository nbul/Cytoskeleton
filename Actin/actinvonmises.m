%% Load distribution data

clear y_forfit1 y_forfit x y bin counter akappa aSD amu p alpha mu_degrees

%% Assign memory
stretch_factor = 0.5;
y_forfit1 = a_added_norm;
counter = size(y_forfit1);
y_forfit = zeros(1,counter(2)-1);
akappa = zeros(1,counter(2)-1);
aSD = zeros(1,counter(2)-1);
amu = zeros(1,counter(2)-1);
x = y_forfit1(:,1);
x = x /90 * pi;
bin = 2*(pi/length(x));

%% Von Mises Fit
for i=2:counter(2)
y_forfit = y_forfit1(:,i);
y = y_forfit;
akappa(i-1) = circ_kappa(x,y,bin);
aSD(i-1) = sqrt(1/akappa(i-1)) * (180/pi) * stretch_factor;
amu(i-1) = circ_mean(x,y);
[p, alpha] = circ_vmpdf(x,amu(i-1),akappa(i-1),x);
amu(i-1) = amu(i-1) * (180/pi) * stretch_factor;
end
