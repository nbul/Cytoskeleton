%% Load distribution data

clear y_forfit1 y_forfit x y bin counter kappa SD mu p alpha mu_degrees

%% Assign memory
stretch_factor = 0.5;
y_forfit1 = m_added_norm;
counter = size(y_forfit1);
y_forfit = zeros(1,counter(2)-1);
kappa = zeros(1,counter(2)-1);
SD = zeros(1,counter(2)-1);
mu = zeros(1,counter(2)-1);
x = y_forfit1(:,1);
x = x /90 * pi;
bin = 0.14;

%% Von Mises Fit
for i=2:counter(2)
y_forfit = y_forfit1(:,i);
y = y_forfit;
kappa(i-1) = circ_kappa(x,y,bin);
SD(i-1) = sqrt(1/kappa(i-1)) * (180/pi) * stretch_factor;
mu(i-1) = circ_mean(x,y);
[p, alpha] = circ_vmpdf(x,mu(i-1),kappa(i-1));
mu(i-1) = mu(i-1) * (180/pi) * stretch_factor;
end






    








