% Secondary leakage inductance with end turns.
%   1  2D harmonic solve, linear core: leakage = L2 - N2*M.
%   2  Energy outside the iron in the 2D model and in an axisymmetric model of the same coil on a
%      round limb, K = L/(mu0 N^2 MLT).
%   3  The 2D model only sees the straight sides, so scale to the real mean turn with the K ratio
%      and compare with L = mu0 N^2 MLT (a + b/3)/h.
% Uses models/ct_model_ac.fem, ct_leak_2d.fem and ct_leak_axi.fem.

1;

function W = outside_iron_energy(groups)
  mo_clearblock;
  for g = groups
    mo_groupselectblock(g);
  end
  W = mo_blockintegral(2);
  mo_clearblock;
end

ct_setup;
ct_params;

openfemm(1);

opendocument(fullfile(modeldir, 'ct_model_ac.fem'));
mi_modifycircprop('Primary', 1, 100);   mi_modifycircprop('Secondary', 1, 0);
mi_analyze(1);  mi_loadsolution;
cp = mo_getcircuitproperties('Secondary');  M = abs(cp(3))/100;
mo_close;
mi_modifycircprop('Primary', 1, 0);     mi_modifycircprop('Secondary', 1, 1);
mi_analyze(1);  mi_loadsolution;
cp = mo_getcircuitproperties('Secondary');  L2 = abs(cp(3));
mo_close;
mi_close;
Lm = N2*M;  Ls2d = L2 - Lm;
fprintf('Part 1 (2D, harmonic solver, linear core)\n');
fprintf('  M = %.3e H, L2 = %.3f mH, Lm = N2*M = %.3f mH, leakage L2 - Lm = %.3f mH (%.1f %% of Lm)\n', ...
        M, L2*1e3, Lm*1e3, Ls2d*1e3, 100*Ls2d/Lm);

opendocument(fullfile(modeldir, 'ct_leak_2d.fem'));
mi_modifycircprop('Primary', 1, 0);  mi_modifycircprop('Secondary', 1, 1);
mi_analyze(1);  mi_loadsolution;
W2d = outside_iron_energy([3 4 5]);
mo_close;
mi_close;
L2d = 2*W2d;
mlt_2d = 2*t;
K2d = (L2d/2)/(mu0*N_limb^2*mlt_2d*1e-3);
fprintf('Part 2 (leakage energy, mu_r = 1e5, no gap)\n');
fprintf('  2D: energy outside iron %.3e J -> %.3f mH for the two coils, K = %.4f\n', W2d, L2d*1e3, K2d);

rl = sqrt(c*t/pi);
opendocument(fullfile(modeldir, 'ct_leak_axi.fem'));
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
