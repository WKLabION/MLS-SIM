function flux = inv_saturate(emmition_rate, lifetime, pulse_width)
% FLUX = saturate(EMMITION_RATE, LIFETIME, PULSE_WIDTH)
%
% solves the emmition power value given specific excitation power with
% flurophore lifetime and pulse width

if nargin == 1
    lifetime = 3.5; % nano seconds
    pulse_width = 0.64; % nano seconds
end

f = @(x) saturate(x, lifetime, pulse_width);

% flux = fzero(@(x)(f(x)-emmition_rate), [0, 30]);

xx = linspace(0, 30, 1000);
yy = f(xx);
flux = interp1(yy, xx, emmition_rate, 'linear', 35);


end