%% Load distribution data

clear y_forfit1 y_forfit x y bin counter kappa SD mu p alpha mu_degrees

%% Assign memory
[file,path] = uigetfile({'*.csv','All CSV files'});
file2 = csvread([path,file]);
stretch_factor = 0.5;
y_forfit1 = file2;
counter = size(y_forfit1);
y_forfit = zeros(1,counter(2)-1);
kappa = zeros(1,counter(2)-1);
SD = zeros(1,counter(2)-1);
mu = zeros(1,counter(2)-1);
x = y_forfit1(:,1);
x = x /90 * pi;
p = zeros(length(x),counter(2)-1);
p2 = zeros(length(x),counter(2)-1);
bin = 2*(pi/length(x));
ecc = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.95, 0.98, 0.98];
b = 1./sqrt(1-ecc .* ecc);
%% Von Mises Fit
for i=2:counter(2)
y_forfit = y_forfit1(:,i);
y = y_forfit;
kappa(i-1) = circ_kappa(x,y,bin);
SD(i-1) = sqrt(1/kappa(i-1)) * (180/pi) * stretch_factor;
mu(i-1) = circ_mean(x,y);
[p(:,i-1), alpha] = circ_vmpdf(x,mu(i-1),kappa(i-1),x);

p2(:,i-1) = 1./(sin((x-mu(i-1))* stretch_factor + pi/2).*sin((x-mu(i-1))* stretch_factor + pi/2)...
    + b(i-1)^2*cos((x-mu(i-1)) * stretch_factor + pi/2).*cos((x-mu(i-1))* stretch_factor + pi/2));
p2(:,i-1) = p2(:,i-1)/sum(p2(:,i-1));
mu(i-1) = mu(i-1) * (180/pi) * stretch_factor;
end


