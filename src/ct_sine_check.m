% Re-solves FEMM at points taken from the loaded runs of ct_load_sim.m and compares flux linkage
% and limb B with the interpolated table. Set I1rms for the 500 A case.

1;

function r = solve_point(i1, i2, xa, xb, cw)
  mi_modifycircprop('Primary', 1, i1);
  mi_modifycircprop('Secondary', 1, i2);
  mi_analyze(1);
  mi_loadsolution;
  cp = mo_getcircuitproperties('Secondary');
  pa = mo_getpointvalues(xa, 0);
  pb = mo_getpointvalues(xb, 0);
  mo_close;
  r = [cp(3), (pa(1) - pb(1))/(cw*1e-3)];
end

ct_setup;
ct_params;
if ~exist('I1rms', 'var'), I1rms = I1_design; end
load(sprintf('ct_check_pts_%dA.dat', I1rms));
xa = bw/2;  xb = bw/2 + c;

openfemm(1);
ct_draw_model(0, 0, 0, []);
mi_saveas(fullfile(modeldir, 'ct_check.fem'));

fprintf('RL [ohm]   i1 [A]   i2 [A]   lambda table   lambda FEMM   error [%%]   B table   B FEMM [T]\n');
err = zeros(size(chk_i, 1), 1);
for k = 1:size(chk_i, 1)
  r = solve_point(chk_i(k,1), chk_i(k,2), xa, xb, c);
  err(k) = 100*abs(r(1) - chk_lam(k))/max(abs(chk_lam));
  fprintf('%7.2f %9.1f %8.3f %12.5f %13.5f %10.3f %9.3f %9.3f\n', ...
          chk_R(k), chk_i(k,1), chk_i(k,2), chk_lam(k), r(1), err(k), chk_B(k), r(2));
end
fprintf('worst flux linkage error %.3f %% of peak\n', max(err));

mi_close;
closefemm;
