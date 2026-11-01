% 50 Hz line current stepped through FEMM from Octave (the magnetostatic solver has no time axis).
% The model has no hysteresis, so only a quarter cycle is solved and the rest follows from
% lambda(-i) = -lambda(i) and lambda(T/2 - t) = lambda(t).
%
% Writes ct_oc_<I>A.dat (open circuit cycle, B map at the current peak) and
% ct_table_<I>A.dat (flux linkage over a grid of line and secondary current).
% Set I1rms before running, default 300:  octave-cli --eval "I1rms=500; ct_sine_femm"

1;

function r = solve_point(i1, i2, xa, xb, cw)
  mi_modifycircprop('Primary', 1, i1);
  mi_modifycircprop('Secondary', 1, i2);
  mi_analyze(1);
  mi_loadsolution;
  cp = mo_getcircuitproperties('Secondary');
  pa = mo_getpointvalues(xa, 0);
  pb = mo_getpointvalues(xb, 0);
  r = [cp(3), (pa(1) - pb(1))/(cw*1e-3)];
end

ct_setup;
ct_params;
if ~exist('I1rms', 'var'), I1rms = I1_design; end
I1pk = I1rms*sqrt(2);
T = 1/f;
wx = bw/2;  wy = hw/2;  ox = wx + c;  oy = wy + c;
xa = wx;  xb = ox;

Nt = 48;
if I1rms <= 300, n1 = 13; n2 = 13; else, n1 = 15; n2 = 13; end
I2max = 1.25*I1pk/N2;

openfemm(1);
ct_draw_model(0, 0, 0, []);
mi_saveas(fullfile(modeldir, 'ct_model_dc.fem'));

nq = Nt/4;
iq = I1pk*sin(2*pi*(0:nq)'/Nt);
lam_q = zeros(nq+1, 1);  B_q = zeros(nq+1, 1);
fprintf('open circuit %g A rms, %d solves\n', I1rms, nq+1);  fflush(stdout);
tic;
for k = 0:nq
  r = solve_point(iq(k+1), 0, xa, xb, c);
  lam_q(k+1) = r(1);  B_q(k+1) = r(2);
  if k == nq
    [X, Y] = meshgrid(-ox+0.5:1:ox-0.5, -oy+0.5:1:oy-0.5);
    ok = ~(abs(X) < wx & abs(Y) < wy) & ~(abs(X) < gap/2 + 0.1 & abs(Y) > wy);
    xs = X(ok);  ys = Y(ok);  Bm = zeros(size(xs));
    for m = 1:numel(xs)
      pv = mo_getpointvalues(xs(m), ys(m));
      Bm(m) = hypot(pv(2), pv(3));
    end
    Bmap = [xs ys Bm];
  end
  if k == 2, fprintf('  %.1f s per solve\n', toc/3); fflush(stdout); end
  if k < nq, mo_close; end
end
mo_close;

time = (0:Nt-1)'*T/Nt;
i1 = I1pk*sin(2*pi*(0:Nt-1)'/Nt);
lam_oc = zeros(Nt, 1);  B_oc = zeros(Nt, 1);
for m = 0:Nt-1
  if m <= nq,          idx = m;          sg = 1;
  elseif m <= 2*nq,    idx = 2*nq - m;   sg = 1;
  elseif m <= 3*nq,    idx = m - 2*nq;   sg = -1;
  else,                idx = Nt - m;     sg = -1;
  end
  lam_oc(m+1) = sg*lam_q(idx+1);
  B_oc(m+1) = sg*B_q(idx+1);
end
save('-text', sprintf('ct_oc_%dA.dat', I1rms), 'f', 'I1pk', 'Nt', 'time', 'i1', 'lam_oc', 'B_oc', 'Bmap');
fprintf('  done in %.0f s\n', toc);  fflush(stdout);

% the model is odd in (i1, i2), so i1 >= 0 is enough
i1g = linspace(0, I1pk, n1);
i2g = linspace(-I2max, I2max, n2);
Lam = zeros(n1, n2);  Bav = zeros(n1, n2);
fprintf('table %d x %d solves\n', n1, n2);  fflush(stdout);
for a = 1:n1
  for b = 1:n2
    r = solve_point(i1g(a), i2g(b), xa, xb, c);
    mo_close;
    Lam(a,b) = r(1);  Bav(a,b) = r(2);
  end
  fprintf('  row %d of %d done\n', a, n1);  fflush(stdout);
end
save('-text', sprintf('ct_table_%dA.dat', I1rms), 'f', 'I1pk', 'i1g', 'i2g', 'Lam', 'Bav');

mi_close;
closefemm;
fprintf('done\n');
