% Hand checks: skin and proximity effect, stray field of the other phases, mass budget.

clear; clc;
ct_setup;
ct_params;

% Dowell, round wire in 3 layers
Tw = 75;
rho_cu = rho_cu20*(1 + alpha_cu*(Tw - 20));
delta_cu = sqrt(2*rho_cu/(w*mu0));
x = d_bare*1e-3/delta_cu;
Delta = (pi/4)^0.75*x*sqrt(d_bare/od);
m = 3;
Fr = 1 + (5*m^2 - 1)/45*Delta^4;
fprintf('Winding at %g C: skin depth %.1f mm, wire %.3f mm (d/delta = %.3f), Dowell R_ac/R_dc for %d layers = %.6f\n', ...
        Tw, delta_cu*1e3, d_bare, x, m, Fr);

rho_al = 2.83e-8*(1 + 0.0040*(75 - 20));
delta_al = sqrt(2*rho_al/(w*mu0));
fprintf('ACSR aluminium at 75 C: skin depth %.1f mm, strand diameter 4.72 mm, ring thickness %.2f mm\n', ...
        delta_al*1e3, (dcond - dsteel)/2);

% both neighbours assumed to add, which is an upper bound
Ipk = I1_design*sqrt(2);
for d_ph = [3 4.5]
  Bext = 2*(mu0*Ipk/(2*pi*d_ph));
  wx = bw/2;  ox = wx + c;
  A_core = 2*ox*1e-3*t*1e-3;
  Phi = Bext*A_core;
  Eext = w*N2*Phi/sqrt(2);
  fprintf('Phase spacing %.1f m: external field <= %.1f uT, induced voltage <= %.1f mV rms (against 16 V), core flux %.2e Wb against %.2e Wb\n', ...
          d_ph, Bext*1e6, Eext*1e3, Phi, 1.4*c*t*1e-6);
end

wy = hw/2;  oy = wy + c;  wx = bw/2;  ox = wx + c;
a_core = 4*ox*oy - 4*wx*wy - 2*gap*c;
m_steel = rho_fe*stack_f*a_core*t*1e-9;
tube = ((c+4)*(t+4) - c*t)*hw;
flange_v = 2*(((c+8)*(t+8) - c*t)*1.0);
m_bob = 2*(tube + flange_v)*1.24e-3*1e-3;     % PLA, 1.24 g/cm3
fprintf('\nMass: steel %.1f g (stacking factor %.3f), copper %.1f g, bobbins %.1f g, total %.1f g\n', ...
        m_steel*1e3, stack_f, m_cu*1e3, m_bob*1e3, (m_steel + m_cu + m_bob)*1e3);
fprintf('Core outer %g x %g mm, window %g x %g mm, limb %g mm, stack %g mm; clear width at full coil height %.2f mm\n', ...
        2*ox, 2*oy, bw, hw, c, t, bw - 2*(wall + 2*od));
