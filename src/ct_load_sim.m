% Loaded runs of the CT (Octave only, FEMM is not called). Uses ct_oc_<I>A.dat and ct_table_<I>A.dat
% from ct_sine_femm.m.
%
% The secondary is stepped in time with
%     d(lambda2)/dt = -R1*i2 + vt,    lambda2 = F(i1, i2) from the FEMM table
% i2 is the FEMM circuit current, so the current leaving the winding is -i2.
%     resistor:  vt = -RL*i2
%     bridge:    vt = +-(Vdc + 2*Vf + 2*rd*|i2|) while two diodes conduct, i2 = 0 otherwise
% Core loss is a resistance Rfe across the induced voltage, so F is evaluated at i1 - N2*e/Rfe.
%
% Set I1rms (default 300) and Tw, the winding temperature for R1 (default 75 C), before running.
% Writes results_<I>A.txt, ct_res_<I>A.dat, ct_check_pts_<I>A.dat and the figures.

1;

function d = time_derivative(x, w)
  n = numel(x);
  X = fft(x(:));
  k = [0:n/2-1, 0, -(n/2-1):-1]';
  d = real(ifft(1i*w*k.*X));
end

function v = table_value(G, i1, i2, tb)
  u = (i1 - tb.i1min)/tb.d1;
  s = (i2 - tb.i2min)/tb.d2;
  u = min(max(u, 0), size(G,1) - 1 - 1e-9);
  s = min(max(s, 0), size(G,2) - 1 - 1e-9);
  a = floor(u);  b = floor(s);
  fu = u - a;    fs = s - b;
  v = (1-fu)*((1-fs)*G(a+1,b+1) + fs*G(a+1,b+2)) + fu*((1-fs)*G(a+2,b+1) + fs*G(a+2,b+2));
end

function r = simulate_load(RL, R1, Rfe, tb, w, I1pk, ncycles, spc)
  dt = (2*pi/w)/spc;  n = ncycles*spc;
  Rt = R1 + RL;  kap = tb.N2*Rt/Rfe;
  rhs = @(t, y) (-Rt*y - table_value(tb.G1, I1pk*sin(w*t) + kap*y, y, tb)*I1pk*w*cos(w*t)) ...
                / (table_value(tb.G2, I1pk*sin(w*t) + kap*y, y, tb) + kap*table_value(tb.G1, I1pk*sin(w*t) + kap*y, y, tb));
  x = zeros(n+1, 1);
  for k = 1:n
    t0 = (k-1)*dt;
    k1 = rhs(t0, x(k));  k2 = rhs(t0+dt/2, x(k)+dt/2*k1);
    k3 = rhs(t0+dt/2, x(k)+dt/2*k2);  k4 = rhs(t0+dt, x(k)+dt*k3);
    x(k+1) = x(k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
  end
  last = (n-spc+2):(n+1);  prev = (n-2*spc+2):(n-spc+1);
  r = finish_run(x(last), -RL*x(last), max(abs(x(last) - x(prev))), R1, Rfe, tb, w, I1pk, n, spc, dt);
  r.iL = -r.i2;
  r.P = mean(RL*r.i2.^2);
  r.V = sqrt(mean(r.vt.^2));  r.I = sqrt(mean(r.iL.^2));
end

function r = simulate_bridge(Vdc, R1, Rfe, Vf, rd, tb, w, I1pk, ncycles, spc)
  dt = (2*pi/w)/spc;  n = ncycles*spc;
  Vth = Vdc + 2*Vf;
  x = zeros(n+1, 1);  vt = zeros(n+1, 1);  mode = 0;
  [vt(1), mode] = open_circuit_voltage(0, I1pk, w, Rfe, tb, Vth);
  for k = 1:n
    t0 = (k-1)*dt;
    if mode ~= 0
      s = mode;
      rhs = @(t, y) bridge_rhs(t, y, s, R1, rd, Vth, Rfe, tb, w, I1pk);
      k1 = rhs(t0, x(k));  k2 = rhs(t0+dt/2, x(k)+dt/2*k1);
      k3 = rhs(t0+dt/2, x(k)+dt/2*k2);  k4 = rhs(t0+dt, x(k)+dt*k3);
      x(k+1) = x(k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
      if s*(-x(k+1)) <= 0, x(k+1) = 0;  mode = 0; end
    end
    if mode == 0
      [vt(k+1), mode] = open_circuit_voltage(k*dt, I1pk, w, Rfe, tb, Vth);
    else
      vt(k+1) = mode*(Vth + 2*rd*abs(x(k+1)));
    end
  end
  last = (n-spc+2):(n+1);  prev = (n-2*spc+2):(n-spc+1);
  r = finish_run(x(last), vt(last), max(abs(x(last) - x(prev))), R1, Rfe, tb, w, I1pk, n, spc, dt);
  r.Idc = mean(abs(r.i2));
  r.Pdc = Vdc*r.Idc;
  r.Pdiode = 2*Vf*r.Idc + 2*rd*mean(r.i2.^2);
  r.Iac = sqrt(mean(r.i2.^2));  r.Ipk = max(abs(r.i2));
  r.Vac = sqrt(mean(r.vt.^2));
  r.cond = mean(r.i2 ~= 0);
end

function [v, mode] = open_circuit_voltage(t, I1pk, w, Rfe, tb, Vth)
  di1 = I1pk*w*cos(w*t);
  e0 = table_value(tb.G1, I1pk*sin(w*t), 0, tb)*di1;
  v = table_value(tb.G1, I1pk*sin(w*t) - tb.N2*e0/Rfe, 0, tb)*di1;
  mode = 0;
  if v > Vth, mode = 1; elseif v < -Vth, mode = -1; end
end

function d = bridge_rhs(t, y, s, R1, rd, Vth, Rfe, tb, w, I1pk)
  e = -R1*y + s*(Vth + 2*rd*abs(y));
  i1e = I1pk*sin(w*t) - tb.N2*e/Rfe;
  F1 = table_value(tb.G1, i1e, y, tb);
  F2 = table_value(tb.G2, i1e, y, tb);
  d = (e - F1*I1pk*w*cos(w*t))/(F2 + F1*tb.N2*(R1 + 2*rd)/Rfe);
end

function r = finish_run(x, vt, settle, R1, Rfe, tb, w, I1pk, n, spc, dt)
  r.t = ((0:spc-1)*dt)';
  r.i2 = x;  r.vt = vt;  r.settle = settle;
  r.i1 = I1pk*sin(w*(n-spc+(1:spc)')*dt);
  r.e = -R1*x + vt;
  r.i1e = r.i1 - tb.N2*r.e/Rfe;
  r.lam = zeros(spc, 1);  r.B = zeros(spc, 1);
  for k = 1:spc
    r.lam(k) = table_value(tb.L, r.i1e(k), x(k), tb);
    r.B(k)   = table_value(tb.B, r.i1e(k), x(k), tb);
  end
  r.Erms = sqrt(mean(r.e.^2));
  r.out_of_table = max(abs(x)) > tb.i2max;
end

function P = steel_loss(Bwave, Bmap, B_ref, stack_f, m_steel, f)
  % local peak B follows the open-circuit map, scaled to the waveform peak
  shape = Bwave/stack_f;
  Bloc = Bmap(:,3)/stack_f*(max(abs(Bwave))/B_ref);
  Bloc = max(Bloc, 0.05);
  p = ct_core_loss(Bloc, shape, f);
  P = m_steel*mean(p);
end

% ------------------------------------------------------------------ main
if ~exist('I1rms', 'var'), I1rms = 300; end
if ~exist('Tw', 'var'), Tw = 75; end
ct_setup;
ct_params;
diary off;
resfile = sprintf('results_%dA.txt', I1rms);
if exist(resfile, 'file'), delete(resfile); end
diary(resfile);

load(sprintf('ct_oc_%dA.dat', I1rms));
load(sprintf('ct_table_%dA.dat', I1rms));
R1 = R1_20*(1 + alpha_cu*(Tw - 20));
wx = bw/2;  wy = hw/2;  ox = wx + c;  oy = wy + c;
a_core = 4*ox*oy - 4*wx*wy - 2*gap*c;
m_steel = rho_fe*stack_f*a_core*t*1e-9;
fprintf('Line current %g A rms (%.1f A peak), %g Hz, N2 = %d\n', I1pk/sqrt(2), I1pk, f, N2);
fprintf('R1 = %.3f ohm at %g C (%.3f at 20 C), steel mass %.1f g\n', R1, Tw, R1_20, m_steel*1e3);

% ---------- open circuit ----------
emf = time_derivative(lam_oc, w);
E = fft(emf)/Nt*2;  Lf = fft(lam_oc)/Nt*2;
Voc = abs(E(2))/sqrt(2);
thd = sqrt(sum(abs(E(3:Nt/2)).^2))/abs(E(2))*100;
Lm = N2*abs(Lf(2))/I1pk;
Bpk = max(abs(B_oc));
Bsteel = Bpk/stack_f;
Hs = interp1(bh(:,1), bh(:,2), Bsteel, 'linear');
fprintf('\nOpen circuit\n');
fprintf('  lambda2 peak %.4f Wb, limb B peak %.3f T (steel %.3f T, H = %.0f A/m, mu_r = %.0f)\n', ...
        max(abs(lam_oc)), Bpk, Bsteel, Hs, Bsteel/(mu0*Hs));
fprintf('  Voc %.2f V rms (%.2f V peak), distortion %.1f %%, Lm %.2f mH\n', Voc, max(abs(emf)), thd, Lm*1e3);
Bmax_local = max(Bmap(:,3))/stack_f;
fprintf('  highest local B in the core at the current peak %.3f T (steel)\n', Bmax_local);
Pfe_oc = steel_loss(B_oc, Bmap, Bpk, stack_f, m_steel, f);
Pfe_oc_avg = m_steel*ct_core_loss(Bsteel, B_oc, f);
[~, Ph1, Pe1] = ct_core_loss(Bsteel, B_oc, f);
fprintf('  core loss %.2f W (local B map), %.2f W (limb B everywhere); hysteresis %.2f, eddy %.2f W/kg\n', ...
        Pfe_oc, Pfe_oc_avg, Ph1, Pe1);
Rfe = Voc^2/Pfe_oc;
fprintf('  equivalent core loss resistance Rfe = %.1f ohm (secondary side)\n', Rfe);

% ---------- (i1, i2) table: mirror to negative i1, then resample ----------
i1f_ = [-fliplr(i1g(2:end)), i1g];
Lf2 = [-Lam(end:-1:2, end:-1:1); Lam];
Bf2 = [-Bav(end:-1:2, end:-1:1); Bav];
Nd = 241;
i1d = linspace(-I1pk, I1pk, Nd)';
i2d = linspace(i2g(1), i2g(end), Nd);
tb.L = interp2(i2g, i1f_, Lf2, i2d, i1d, 'spline');
tb.B = interp2(i2g, i1f_, Bf2, i2d, i1d, 'spline');
tb.d1 = i1d(2) - i1d(1);  tb.d2 = i2d(2) - i2d(1);
tb.i1min = i1d(1);  tb.i2min = i2d(1);  tb.i2max = i2d(end);  tb.N2 = N2;
[tb.G2, tb.G1] = gradient(tb.L, tb.d2, tb.d1);
lam_tab = zeros(Nt, 1);
for k = 1:Nt, lam_tab(k) = table_value(tb.L, i1(k), 0, tb); end
fprintf('\nTable %d x %d FEMM solves; interpolated table vs open circuit: worst error %.3f %% of peak\n', ...
        numel(i1g), numel(i2g), 100*max(abs(lam_tab - lam_oc))/max(abs(lam_oc)));

ncyc = 10;  spc = 1000;

% ---------- resistor at the design burden ----------
ld = simulate_load(RL_report, R1, Rfe, tb, w, I1pk, ncyc, spc);
Pfe_ld = steel_loss(ld.B, Bmap, Bpk, stack_f, m_steel, f);
Pcu_ld = R1*mean(ld.i2.^2);
fprintf('\nResistor RL = %.2f ohm\n', RL_report);
fprintf('  load %.2f V, %.3f A, %.2f W;  winding loss %.2f W, core loss %.2f W;  limb B peak %.3f T\n', ...
        ld.V, ld.I, ld.P, Pcu_ld, Pfe_ld, max(abs(ld.B)));
fprintf('  cycle repeats to %.1e A\n', ld.settle);

% ---------- burden sweep ----------
Rs = [1 2 4 6 8 9.88 12 16 24 40];
sw = zeros(numel(Rs), 5);
for k = 1:numel(Rs)
  r = simulate_load(Rs(k), R1, Rfe, tb, w, I1pk, 8, 500);
  sw(k,:) = [r.P, r.V, r.I, max(abs(r.B)), R1*mean(r.i2.^2)];
end
[~, kb] = max(sw(:,1));
lo = Rs(max(kb-1, 1));  hi = Rs(min(kb+1, numel(Rs)));
Ropt = fminbnd(@(x) -simulate_load(x, R1, Rfe, tb, w, I1pk, 8, 500).P, lo, hi, optimset('TolX', 0.02));
opt = simulate_load(Ropt, R1, Rfe, tb, w, I1pk, ncyc, spc);
Pfe_opt = steel_loss(opt.B, Bmap, Bpk, stack_f, m_steel, f);
Pcu_opt = R1*mean(opt.i2.^2);
fprintf('\nBurden sweep\n   RL [ohm]    P [W]   V [V]   I [A]   Bpk [T]   Pcu [W]\n');
for k = 1:numel(Rs), fprintf('   %6.2f   %7.2f %7.2f %7.3f %8.3f %8.2f\n', Rs(k), sw(k,:)); end
fprintf('Best burden %.2f ohm: P = %.2f W, V = %.2f V, I = %.3f A, B peak %.3f T, winding %.2f W, core %.2f W\n', ...
        Ropt, opt.P, opt.V, opt.I, max(abs(opt.B)), Pcu_opt, Pfe_opt);
fprintf('Efficiency at the best burden: %.1f %% of the power taken from the line (load / (load + winding + core))\n', ...
        100*opt.P/(opt.P + Pcu_opt + Pfe_opt));

% ---------- bridge rectifier and dc link ----------
Vs = [4 6 8 10 12 14 16 18 20 24];
Pdc = zeros(size(Vs));  Idc = Pdc;  Iac = Pdc;  Ipk = Pdc;  Pd = Pdc;  Pcb = Pdc;  cnd = Pdc;  Vac = Pdc;
for k = 1:numel(Vs)
  r = simulate_bridge(Vs(k), R1, Rfe, Vf, rd, tb, w, I1pk, 8, 1000);
  Pdc(k) = r.Pdc;  Idc(k) = r.Idc;  Iac(k) = r.Iac;  Ipk(k) = r.Ipk;
  Pd(k) = r.Pdiode;  Pcb(k) = R1*mean(r.i2.^2);  cnd(k) = r.cond;  Vac(k) = r.Vac;
end
[~, kbr] = max(Pdc);
Vbest = fminbnd(@(v) -simulate_bridge(v, R1, Rfe, Vf, rd, tb, w, I1pk, 8, 1000).Pdc, ...
                Vs(max(kbr-1,1)), Vs(min(kbr+1,numel(Vs))), optimset('TolX', 0.1));
brb = simulate_bridge(Vbest, R1, Rfe, Vf, rd, tb, w, I1pk, 10, 2000);
Pfe_br = steel_loss(brb.B, Bmap, Bpk, stack_f, m_steel, f);
fprintf('\nBridge rectifier (Schottky Vf = %.2f V, rd = %.2f ohm) into a dc link held at Vdc\n', Vf, rd);
fprintf('   Vdc [V]  Pdc [W]  Idc [A]  Iac rms  Iac peak  conduction  winding W  diode W\n');
for k = 1:numel(Vs)
  fprintf('   %6.1f %8.2f %8.3f %8.3f %8.3f %9.0f %% %9.2f %8.2f\n', Vs(k), Pdc(k), Idc(k), Iac(k), Ipk(k), 100*cnd(k), Pcb(k), Pd(k));
end
fprintf('Best dc link voltage %.1f V: Pdc = %.2f W, Idc = %.3f A, ac current %.3f A rms / %.3f A peak, conduction %.0f %% of the cycle\n', ...
        Vbest, brb.Pdc, brb.Idc, brb.Iac, brb.Ipk, 100*brb.cond);
fprintf('   winding loss %.2f W, diode loss %.2f W, core loss %.2f W, equivalent ac resistance %.2f ohm\n', ...
        R1*mean(brb.i2.^2), brb.Pdiode, Pfe_br, brb.Vac/brb.Iac);

% ---------- charger and battery ----------
fprintf('\nCharger: buck stage %.0f %% efficient, 2S pack, Rbat = %.2f ohm, input voltage regulated\n', 100*eta_buck, Rbat);
fprintf('   Vbat   Vdc used   Pdc [W]   Pbat [W]   Ichg [A]   charge time for %.1f Ah [h]\n', Cbat);
chg = zeros(numel(Vbat), 5);
for q = 1:numel(Vbat)
  Vmin = Vbat(q) + 1.5;                                  % buck headroom
  Vu = max(Vbest, Vmin);
  rr = simulate_bridge(Vu, R1, Rfe, Vf, rd, tb, w, I1pk, 8, 1000);
  Pb = eta_buck*rr.Pdc;
  Ic = (-Vbat(q) + sqrt(Vbat(q)^2 + 4*Rbat*Pb))/(2*Rbat);
  th = Cbat*mean(Vbat)/Pb;
  chg(q,:) = [Vu, rr.Pdc, Pb, Ic, th];
  fprintf('   %4.1f V  %7.1f V  %8.2f  %9.2f  %9.2f   %8.1f\n', Vbat(q), Vu, rr.Pdc, Pb, Ic, th);
end

% ---------- saved results ----------
res.I1rms = I1rms;  res.Tw = Tw;  res.R1 = R1;  res.Voc = Voc;  res.thd = thd;  res.Lm = Lm;
res.Bpk_oc = Bpk;  res.Bsteel_oc = Bsteel;  res.Bmax_local = Bmax_local;
res.Pfe_oc = Pfe_oc;  res.Pfe_oc_avg = Pfe_oc_avg;  res.Rfe = Rfe;  res.m_steel = m_steel;
res.RL = RL_report;  res.ld = [ld.V ld.I ld.P Pcu_ld Pfe_ld max(abs(ld.B))];
res.Ropt = Ropt;  res.opt = [opt.V opt.I opt.P Pcu_opt Pfe_opt max(abs(opt.B))];
res.sw = [Rs(:) sw];
res.Vs = Vs;  res.bridge = [Vs(:) Pdc(:) Idc(:) Iac(:) Ipk(:) cnd(:) Pcb(:) Pd(:)];
res.Vbest = Vbest;  res.brb = [brb.Pdc brb.Idc brb.Iac brb.Ipk brb.cond R1*mean(brb.i2.^2) brb.Pdiode Pfe_br brb.Vac/brb.Iac];
res.chg = chg;  res.Vbat = Vbat;
save('-text', sprintf('ct_res_%dA.dat', I1rms), 'res');

pick = 1:spc/8:spc;
chk_i = [ld.i1(pick), ld.i2(pick); opt.i1(pick), opt.i2(pick)];
chk_lam = [ld.lam(pick); opt.lam(pick)];
chk_B = [ld.B(pick); opt.B(pick)];
chk_R = [RL_report*ones(numel(pick),1); Ropt*ones(numel(pick),1)];
save('-text', sprintf('ct_check_pts_%dA.dat', I1rms), 'chk_i', 'chk_lam', 'chk_B', 'chk_R');

% ---------- figures ----------
graphics_toolkit('gnuplot');
tm = time*1e3;  tl = ld.t*1e3;
fig = figure('visible', 'off', 'position', [100 100 1300 750]);
graphics_toolkit(fig, 'gnuplot');
subplot(2,3,1);  plot(tm, i1, 'o-', 'markersize', 3);  grid on;
title('Line current i_1(t)');  xlabel('t [ms]');  ylabel('A');
subplot(2,3,2);  plot(tm, lam_oc, 'o-', 'markersize', 3);  grid on;
title('Open circuit flux linkage');  xlabel('t [ms]');  ylabel('Wb-turn');
subplot(2,3,3);  plot(tm, emf, 'o-', 'markersize', 3);  grid on;
title(sprintf('Open circuit EMF, %.2f V rms', Voc));  xlabel('t [ms]');  ylabel('V');
subplot(2,3,4);  plot(tl, ld.vt, 'b', tl, ld.iL*RL_report, 'r--');  grid on;
title(sprintf('Load %.2f ohm: voltage (blue), current x R (dash)', RL_report));  xlabel('t [ms]');  ylabel('V');
subplot(2,3,5);  plot(tl, RL_report*ld.i2.^2, 'b', tl, ld.P*ones(size(tl)), 'k--');  grid on;
title(sprintf('Load power, mean %.2f W', ld.P));  xlabel('t [ms]');  ylabel('W');
subplot(2,3,6);  plot(tm, B_oc, 'o-', 'markersize', 3, tl, ld.B, 'r');  grid on;
title('Limb B: open circuit (o), loaded (red)');  xlabel('t [ms]');  ylabel('T');
print(fig, sprintf('fig_sine_%dA.png', I1rms), '-dpng', '-r100');

fig2 = figure('visible', 'off', 'position', [100 100 900 420]);
graphics_toolkit(fig2, 'gnuplot');
subplot(1,2,1);  semilogx(Rs, sw(:,1), 'o-', Ropt, opt.P, 'r*', RL_report, ld.P, 'ks');  grid on;
xlabel('R_L [ohm]');  ylabel('P [W]');  title('Load power against burden');
subplot(1,2,2);  semilogx(Rs, sw(:,2), 'o-', Rs, sw(:,3), 's-');  grid on;
xlabel('R_L [ohm]');  title('Voltage [V rms] (o), current [A rms] (s)');
print(fig2, sprintf('fig_sweep_%dA.png', I1rms), '-dpng', '-r100');

fig4 = figure('visible', 'off', 'position', [100 100 1300 420]);
graphics_toolkit(fig4, 'gnuplot');
tb_ = (0:numel(brb.i2)-1)'*(1/f)/numel(brb.i2)*1e3;
subplot(1,3,1);  plot(tb_, brb.i2, 'b', tb_, brb.i1/N2, 'k--');  grid on;
title(sprintf('Bridge, Vdc = %.1f V: secondary current', Vbest));  xlabel('t [ms]');  ylabel('A');
subplot(1,3,2);  plot(tb_, brb.vt, 'r');  grid on;
title('Voltage at the secondary terminals');  xlabel('t [ms]');  ylabel('V');
subplot(1,3,3);  plot(Vs, Pdc, 'o-', Vbest, brb.Pdc, 'r*');  grid on;
title('Power into the dc link');  xlabel('V_{dc} [V]');  ylabel('W');
print(fig4, sprintf('fig_bridge_%dA.png', I1rms), '-dpng', '-r100');

if I1rms ~= I1_design
  c3 = load('ct_oc_300A.dat');
  e3 = time_derivative(c3.lam_oc, w);
  s = I1pk/c3.I1pk;
  V3 = abs(fft(e3)(2))/c3.Nt*2/sqrt(2);
  fprintf('\nSaturation, open circuit, 300 A -> %g A rms (current x%.3f)\n', I1pk/sqrt(2), s);
  fprintf('  B peak %.3f -> %.3f T (x%.3f, a linear core gives x%.3f)\n', ...
          max(abs(c3.B_oc)), Bpk, Bpk/max(abs(c3.B_oc)), s);
  fprintf('  Voc %.2f -> %.2f V rms (x%.3f)\n', V3, Voc, Voc/V3);
  fig3 = figure('visible', 'off', 'position', [100 100 1300 420]);
  graphics_toolkit(fig3, 'gnuplot');
  subplot(1,3,1);  plot(c3.time*1e3, c3.B_oc*s, 'k--', c3.time*1e3, c3.B_oc, 'b', tm, B_oc, 'r');  grid on;
  title('Limb B: 300 A (blue), 500 A (red)');  xlabel('t [ms]');  ylabel('T');
  subplot(1,3,2);  plot(c3.time*1e3, e3*s, 'k--', c3.time*1e3, e3, 'b', tm, emf, 'r');  grid on;
  title('Open circuit EMF (dash: 300 A x 5/3)');  xlabel('t [ms]');  ylabel('V');
  subplot(1,3,3);  plot(i1g, Bav(:, ceil(end/2)), 'r', i1g, i1g*max(abs(c3.B_oc))/c3.I1pk, 'k--');  grid on;
  title('B against i_1 at i_2 = 0 (dash: 300 A slope)');  xlabel('i_1 [A]');  ylabel('T');
  print(fig3, 'fig_saturation.png', '-dpng', '-r100');
end

diary off;
