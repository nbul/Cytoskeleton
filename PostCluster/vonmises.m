%% Assign memory
counter = size(m_added_norm);
y_forfit = zeros(1,counter(2)-1);
kappa = zeros(1,counter(2)-1);
SD = zeros(1,counter(2)-1);
mu = zeros(1,counter(2)-1);
x = m_added_norm(:,1);
x = x /90 * pi;
bin = 2*pi/length(x);

%% Von Mises Fit
for i=2:counter(2)
kappa(i-1) = circ_kappa(x,m_added_norm(:,i),bin);
SD(i-1) = sqrt(1/kappa(i-1)) * (180/pi) * stretch_factor;
mu(i-1) = circ_mean(x,m_added_norm(:,i));
[p, alpha] = circ_vmpdf(x,mu(i-1),kappa(i-1),bin);
mu(i-1) = mu(i-1) * (180/pi) * stretch_factor;
end






    








