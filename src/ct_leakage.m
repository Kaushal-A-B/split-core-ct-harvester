% Secondary leakage inductance including end turns.
%   1  2D harmonic solve with a linear core: M from the line current, L2 from the secondary
%      alone; leakage = L2 - N2*M.
%   2  Leakage energy outside the iron (mu_r = 1e5, no joints) for the 2D coil and for the same
%      coil on a round limb of equal area in an axisymmetric model. K = L/(mu0 N^2 MLT).
%   3  The 2D model sees only the straight sides, so the leakage is scaled to the real mean turn
%      with the K ratio of part 2 and compared with L = mu0 N^2 MLT (a + b/3)/h.

1;

function W = outside_iron_energy(groups)
  mo_clearblock;
  for g = groups
    mo_groupselectblock(g);
  end
  W = mo_blockintegral(2);
  mo_clearblock;
end

function draw_poly(P)
  n = size(P, 1);
  for k = 1:n, mi_addnode(P(k,1), P(k,2)); end
  for k = 1:n
    j = mod(k, n) + 1;
    mi_addsegment(P(k,1), P(k,2), P(j,1), P(j,2));
  end
end

function label(x, y, mat, auto, msize, circ, turns, grp)
  mi_addblocklabel(x, y);
  mi_selectlabel(x, y);
  mi_setblockprop(mat, auto, msize, circ, 0, grp, turns);
  mi_clearselected;
end

ct_setup;
ct_params;
wx = bw/2;  wy = hw/2;  ox = wx + c;

openfemm(1);

mu_lin = 42345;
ct_draw_model(0, 0, f, mu_lin);
mi_saveas(fullfile(modeldir, 'ct_leak_ac.fem'));
mi_modifycircprop('Primary', 1, 100);   mi_modifycircprop('Secondary', 1, 0);
mi_analyze(1);  mi_loadsolution;
cp = mo_getcircuitproperties('Secondary');  M = abs(cp(3))/100;
mo_close;
mi_modifycircprop('Primary', 1, 0);     mi_modifycircprop('Secondary', 1, 1);
mi_analyze(1);  mi_loadsolution;
cp = mo_getcircuitproperties('Secondary');  L2 = abs(cp(3));
mo_close;
Lm = N2*M;  Ls2d = L2 - Lm;
fprintf('Part 1 (2D, harmonic solver, linear core)\n');
fprintf('  M = %.3e H, L2 = %.3f mH, Lm = N2*M = %.3f mH, leakage L2 - Lm = %.3f mH (%.1f %% of Lm)\n', ...
        M, L2*1e3, Lm*1e3, Ls2d*1e3, 100*Ls2d/Lm);

ct_draw_model(0, 0, 0, 1e5, true);
mi_saveas(fullfile(modeldir, 'ct_leak_2d.fem'));
mi_modifycircprop('Primary', 1, 0);  mi_modifycircprop('Secondary', 1, 1);
mi_analyze(1);  mi_loadsolution;
W2d = outside_iron_energy([3 4 5]);
mo_close;
L2d = 2*W2d;
mlt_2d = 2*t;
K2d = (L2d/2)/(mu0*N_limb^2*mlt_2d*1e-3);
fprintf('Part 2 (leakage energy, mu_r = 1e5, no gap)\n');
fprintf('  2D: energy outside iron %.3e J -> %.3f mH for the two coils, K = %.4f\n', W2d, L2d*1e3, K2d);

rl = sqrt(c*t/pi);
rw = rl + bw/2;
Ro = rw + c;
yh = hw/2 + c;
newdocument(0);
mi_probdef(0, 'millimeters', 'axi', 1e-8, 0, 30);
mi_drawrectangle(0, -yh, Ro, yh);
mi_drawrectangle(rl, -hw/2, rw, hw/2);
t2 = 2*od;  t3 = od;  hp = n_part*od;
xc1 = rl + wall;  xc2 = xc1 + t2;  xd2 = xc2 + t3;
Q = [xc1 -Lw/2; xc2 -Lw/2; xc2 Lw/2-hp; xd2 Lw/2-hp; xd2 Lw/2; xc1 Lw/2];
draw_poly(Q);
mi_addmaterial('Fe', 1e5, 1e5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
mi_addmaterial('Coil_wire', 1, 1, 0, 0, 58, 0, 0, 1, 3, 0, 0, 1, d_bare);
mi_getmaterial('Air');
mi_addcircprop('Coil', 1, 1);
label(rl/2, 0, 'Fe', 1, 0, '', 0, 1);
label(rw - 1, 0, 'Air', 1, 0, '', 0, 2);
label(Ro + 20, 0, 'Air', 1, 0, '', 0, 2);
label(xc1 + od, 0, 'Coil_wire', 1, 0, 'Coil', N_limb, 3);
mi_makeABC(7, 3*Ro, 0, 0, 0);
mi_saveas(fullfile(modeldir, 'ct_leak_axi.fem'));
mi_analyze(1);  mi_loadsolution;
Wax = outside_iron_energy([2 3]);
mo_close;
Lax = 2*Wax;
mlt_ax = 2*pi*(rl + wall + 0.75);
Kax = Lax/(mu0*N_limb^2*mlt_ax*1e-3);
fprintf('  axisymmetric (round limb r = %.1f mm): energy %.3e J -> %.3f mH for one coil, K = %.4f\n', rl, Wax, Lax*1e3, Kax);
mi_close;
closefemm;

mlt = len_limb/N_limb;
scale = mlt/mlt_2d*Kax/K2d;
Ls3d = Ls2d*scale;
b_coil = 3*od;
K_formula = (wall + b_coil/3)/Lw;
L_formula = 2*mu0*N_limb^2*mlt*1e-3*K_formula;
fprintf('Part 3\n');
fprintf('  real mean turn %.1f mm against %.1f mm in the 2D model, K ratio axisymmetric/2D = %.2f\n', mlt, mlt_2d, Kax/K2d);
fprintf('  leakage with end turns: %.3f mH x %.2f = %.3f mH (%.1f %% of Lm)\n', Ls2d*1e3, scale, Ls3d*1e3, 100*Ls3d/Lm);
fprintf('  formula mu0 N^2 MLT (a + b/3)/h for two coils: %.3f mH (K = %.4f)\n', L_formula*1e3, K_formula);
save('-text', 'ct_leakage.dat', 'M', 'L2', 'Lm', 'Ls2d', 'Ls3d', 'scale', 'K2d', 'Kax', 'K_formula', 'L_formula', 'mlt');
