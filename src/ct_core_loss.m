function [P, Ph, Pe] = ct_core_loss(Bpk, shape, f)
% Core loss (W/kg) of 0.27 mm M-4. Bpk is the peak flux density, shape one cycle of the B wave.
% Eddy loss follows the mean square of dB/dt, so a flattened wave pays for its steep edges.
% Coefficients from ct_params.m, extrapolated log-log outside 1 to 1.8 T.

ct_params;

n = numel(shape);
x = shape(:)/max(abs(shape(:)));
X = fft(x);
k = [0:n/2-1, 0, -(n/2-1):-1]';
dx = real(ifft(1i*2*pi*f*k.*X));
ratio = mean(dx.^2)/((2*pi*f)^2/2);

a = exp(interp1(log(B_loss), log(a_h), log(Bpk), 'linear', 'extrap'));
b = exp(interp1(log(B_loss), log(b_e), log(Bpk), 'linear', 'extrap'));
Ph = a*f;
Pe = b*f^2*ratio;
P  = Ph + Pe;
end
