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
for fit=2:counter(2)
y_forfit = y_forfit1(:,fit);
y = y_forfit;
kappa(fit-1) = circ_kappa(x,y,bin);
SD(fit-1) = sqrt(1/kappa(fit-1)) * (180/pi) * stretch_factor;
mu(fit-1) = circ_mean(x,y);
[p, alpha] = circ_vmpdf(x,mu(fit-1),kappa(fit-1));
mu(fit-1) = mu(fit-1) * (180/pi) * stretch_factor;
end






    








