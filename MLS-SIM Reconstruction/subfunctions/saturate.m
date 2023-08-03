function emmition_rate = saturate(flux, lifetime, pulse_width)
% EMMITION_RATE = saturate(FLUX, LIFETIME, PULSE_WIDTH)
%
% solves the emmition power value given specific excitation power with
% flurophore lifetime and pulse width


% excitation power in photons per unit time in one cross section
% excitation_power = 1; % photons per nano seconds
% lifetime of flurophore in SSIM paper is 3.5 ns
% lifetime of a GFP flurophore is 2.6 ns
% lifetime = 2.6; % nano seconds
% pulse width of laser in SSIM paper is 0.64 ns
% pulse_width = 0.64; % nano seconds
if nargin == 1
    lifetime = 3.5; % nano seconds
    pulse_width = 0.64; % nano seconds
end
decay_rate = 1 / lifetime;

% emmition is proportional to 
% 1/[1+1/(lifetime*flux_intensity*crosssection)]*{1-exp[-(flux_intensity*crosssection+1/lifetime)*pulse_width]}

pt1 = flux ./ (flux + decay_rate);
pt2 = 1 - exp(-(flux + decay_rate) .* pulse_width);
% in arbitrary unit
emmition_rate = pt1 .* pt2;

% flux_per_lifetime = flux .* lifetime;
% figure;
% plot(flux_per_lifetime, emmition_rate); hold on;
% plot(flux_per_lifetime, pt1); hold on;
% plot(flux_per_lifetime, pt2); hold off;

end