% Two-node (core, winding) lumped thermal model in still air, using the losses from ct_load_sim.m.
% The hot line conductor warms the window-side surfaces. Needs ct_res_300A.dat and ct_res_500A.dat.

1;
clear; clc;
ct_setup;
ct_params;

Ta = 40;                         % deg C
T_line300 = 75;                  % deg C at 300 A, assumed
h_list = [5 10 20];              % W/m2K
h = 10;
k_pla = 0.13;  cp_pla = 1800;  m_pla = 13.8e-3;
Tpla = 60;                       % deg C, PLA softening
Tcu_max = 155;                   % deg C, class F enamel

wx = bw/2;  wy = hw/2;  ox = wx + c;  oy = wy + c;
A_c_out = (2*(2*ox + 2*oy)*t/2 + 2*(4*ox*oy - 4*wx*wy))*1e-6 - 2*Lw*t*1e-6 - 2*2*Lw*c*1e-6;
A_c_win = (2*(2*wx + 2*wy)*t)*1e-6 - 2*Lw*t*1e-6;
A_coil = 2*((2*(c+t) + 2*pi*(wall + 2*od))*Lw + 2*((c+8)*(t+8) - c*t))*1e-6;
A_w_win = 0.30*A_coil;  A_w_out = A_coil - A_w_win;
G_cw = k_pla*(2*(2*(c+t))*Lw*1e-6)/(wall*1e-3);
C_c = rho_fe*stack_f*(4*ox*oy - 4*wx*wy)*t*1e-9*cp_fe;
C_w = m_cu*385 + m_pla*cp_pla;
fprintf('Areas (cm2): core outside %.0f, core window faces %.0f, coils outside %.0f, coils window side %.0f\n', ...
        1e4*A_c_out, 1e4*A_c_win, 1e4*A_w_out, 1e4*A_w_win);
fprintf('Bobbin conduction %.2f W/K, heat capacity core %.0f J/K, winding + bobbin %.0f J/K\n\n', G_cw, C_c, C_w);

cases = {300, 'bridge'; 300, 'resistor'; 500, 'bridge'; 500, 'resistor'};
out = zeros(size(cases, 1), 8);
for q = 1:size(cases, 1)
  I = cases{q,1};
  d = load(sprintf('ct_res_%dA.dat', I));  res = d.res;
  if strcmp(cases{q,2}, 'bridge')
    Pfe = res.brb(8);  Pcu_ref = res.brb(6);  Pload = res.brb(1);
  else
    Pfe = res.opt(5);  Pcu_ref = res.opt(4);  Pload = res.opt(3);
  end
  T_line = Ta + (T_line300 - Ta)*(I/300)^2;
  for hh = h_list
    G_c_out = hh*A_c_out;  G_c_win = hh*A_c_win;  G_w_out = hh*A_w_out;  G_w_win = hh*A_w_win;
    Twin = (Ta + T_line)/2;
    Tc = Ta;  Tw = Ta;
    for it = 1:50
      Pcu = Pcu_ref*(1 + alpha_cu*(Tw - 20))/(1 + alpha_cu*(res.Tw - 20));
      A = [G_c_out + G_c_win + G_cw, -G_cw; -G_cw, G_w_out + G_w_win + G_cw];
      b = [Pfe + G_c_out*Ta + G_c_win*Twin; Pcu + G_w_out*Ta + G_w_win*Twin];
      T = A\b;  Tc = T(1);  Tw = T(2);
    end
    if hh == h, out(q,:) = [I, Pload, Pfe, Pcu, Tc, Tw, Twin, T_line]; end
    fprintf('%g A %-8s h = %2g: load %.1f W, core loss %.2f W, winding loss %.2f W (R1 at %.0f C) -> core %.0f C, winding %.0f C\n', ...
            I, cases{q,2}, hh, Pload, Pfe, Pcu, Tw, Tc, Tw);
  end
end
fprintf('\nR1 in ct_load_sim.m was taken at %g C; rerun it with Tw set if this differs by more than about 10 K.\n', res.Tw);

q = 3;  Pfe = out(q,3);  Pcu0 = out(q,4);  T_line = out(q,8);  Twin = out(q,7);
dt = 5;  n = 3*3600/dt;  Tc = Ta;  Tw = Ta;  hist = zeros(n, 3);
t60 = NaN;
for k = 1:n
  Pcu = Pcu0*(1 + alpha_cu*(Tw - 20))/(1 + alpha_cu*(out(q,6) - 20));
  dTc = (Pfe - h*A_c_out*(Tc - Ta) - h*A_c_win*(Tc - Twin) - G_cw*(Tc - Tw))/C_c;
  dTw = (Pcu - h*A_w_out*(Tw - Ta) - h*A_w_win*(Tw - Twin) - G_cw*(Tw - Tc))/C_w;
  Tc = Tc + dt*dTc;  Tw = Tw + dt*dTw;
  hist(k,:) = [k*dt/60, Tc, Tw];
  if isnan(t60) && Tw >= Tpla, t60 = k*dt/60; end
end
fprintf('500 A, bridge, warm-up from %g C: winding reaches the PLA limit of %g C after %.0f min\n', Ta, Tpla, t60);
fprintf('  after 30 min: core %.0f C, winding %.0f C;  after 3 h: core %.0f C, winding %.0f C\n', ...
        hist(30*60/dt,2), hist(30*60/dt,3), hist(end,2), hist(end,3));

q = 1;  Pfe = out(q,3);  Pcu0 = out(q,4);  Twin = out(q,7);
Tc = Ta;  Tw = Ta;  t60b = NaN;
for k = 1:n
  Pcu = Pcu0*(1 + alpha_cu*(Tw - 20))/(1 + alpha_cu*(out(q,6) - 20));
  dTc = (Pfe - h*A_c_out*(Tc - Ta) - h*A_c_win*(Tc - Twin) - G_cw*(Tc - Tw))/C_c;
  dTw = (Pcu - h*A_w_out*(Tw - Ta) - h*A_w_win*(Tw - Twin) - G_cw*(Tw - Tc))/C_w;
  Tc = Tc + dt*dTc;  Tw = Tw + dt*dTw;
  hist3(k,:) = [k*dt/60, Tc, Tw];
  if isnan(t60b) && Tw >= Tpla, t60b = k*dt/60; end
end
fprintf('300 A, bridge: after 3 h core %.0f C, winding %.0f C; PLA limit %g C reached after %.0f min\n', hist3(end,2), hist3(end,3), Tpla, t60b);

save('-text', 'ct_thermal.dat', 'out', 'hist', 'hist3', 'Ta', 'Tpla', 'Tcu_max', 't60', 't60b');

graphics_toolkit('gnuplot');
fig = figure('visible', 'off', 'position', [100 100 800 450]);
graphics_toolkit(fig, 'gnuplot');
plot(hist(:,1), hist(:,3), 'r', hist(:,1), hist(:,2), 'b', hist3(:,1), hist3(:,3), 'r--', hist3(:,1), hist3(:,2), 'b--', ...
     [0 180], [Tpla Tpla], 'k:');
grid on;  xlabel('time [min]');  ylabel('temperature [C]');
title('Warm-up, bridge load: winding (red), core (blue), 500 A solid, 300 A dashed');
print(fig, 'fig_thermal.png', '-dpng', '-r100');
