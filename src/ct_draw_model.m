function ct_draw_model(xo, yo, freq, mu_lin, nogap)
% Draws the CT in an open FEMM session (planar, front view).
%   xo, yo   conductor axis position in the window (mm)
%   freq     0 magnetostatic, 50 harmonic
%   mu_lin   [] for the M-4 B-H curve, a number for a linear core
%   nogap    true fills the joints with core steel

ct_params;
if nargin < 5, nogap = false; end

wx = bw/2;  wy = hw/2;
ox = wx + c;  oy = wy + c;
g = gap/2;

newdocument(0);
mi_probdef(freq, 'millimeters', 'planar', 1e-8, t, 30);

mi_drawrectangle(-ox, -oy, ox, oy);
mi_drawrectangle(-wx, -wy, wx, wy);

for s = [1 -1]
  mi_addnode(-g, s*wy);  mi_addnode(g, s*wy);
  mi_addnode(-g, s*oy);  mi_addnode(g, s*oy);
  mi_addsegment(-g, s*wy, -g, s*oy);
  mi_addsegment( g, s*wy,  g, s*oy);
end

ro = dcond/2;  rs = dsteel/2;
mi_addnode(xo-ro, yo);  mi_addnode(xo+ro, yo);
mi_addarc(xo-ro, yo, xo+ro, yo, 180, 2.5);
mi_addarc(xo+ro, yo, xo-ro, yo, 180, 2.5);
mi_addnode(xo-rs, yo);  mi_addnode(xo+rs, yo);
mi_addarc(xo-rs, yo, xo+rs, yo, 180, 2.5);
mi_addarc(xo+rs, yo, xo-rs, yo, 180, 2.5);

t2 = 2*od;  t3 = od;  hp = n_part*od;
for side = [1 -1]
  xa1 = wx - wall - t2;  xa2 = wx - wall;  xb1 = xa1 - t3;
  P = [xa2 -Lw/2; xa1 -Lw/2; xa1 Lw/2-hp; xb1 Lw/2-hp; xb1 Lw/2; xa2 Lw/2];
  draw_poly(P.*[side 1]);
  xc1 = ox + wall;  xc2 = xc1 + t2;  xd2 = xc2 + t3;
  Q = [xc1 -Lw/2; xc2 -Lw/2; xc2 Lw/2-hp; xd2 Lw/2-hp; xd2 Lw/2; xc1 Lw/2];
  draw_poly(Q.*[side 1]);
end

if isempty(mu_lin)
  mi_addmaterial('M4', 1, 1, 0, 0, 0, lam_d, 0, stack_f, 0, 0, 0, 0, 0);
  for k = 1:size(bh, 1)
    mi_addbhpoint('M4', bh(k,1), bh(k,2));
  end
else
  mi_addmaterial('M4', mu_lin, mu_lin, 0, 0, 0, lam_d, 0, stack_f, 0, 0, 0, 0, 0);
end
mi_addmaterial('Coil_wire', 1, 1, 0, 0, 58, 0, 0, 1, 3, 0, 0, 1, d_bare);
mi_getmaterial('Air');
mi_getmaterial('Aluminum, 1100');
if freq == 0
  steel_name = '1010 Steel';
  mi_getmaterial(steel_name);
else
  % eddy currents in the nonlinear ACSR steel strands are not worth the mesh at 50 Hz
  steel_name = 'ACSR_core';
  mi_addmaterial(steel_name, 100, 100, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
end

mi_addcircprop('Primary', I1_design*sqrt(2), 1);
mi_addcircprop('Secondary', 0, 1);

xcore = (wx + ox)/2;
label( xcore, 0, 'M4', 1, 0, '', 0, 1);
label(-xcore, 0, 'M4', 1, 0, '', 0, 1);
if nogap
  label(0,  (wy+oy)/2, 'M4', 1, 0, '', 0, 1);
  label(0, -(wy+oy)/2, 'M4', 1, 0, '', 0, 1);
else
  label(0,  (wy+oy)/2, 'Air', 0, 0.05, '', 0, 2);
  label(0, -(wy+oy)/2, 'Air', 0, 0.05, '', 0, 2);
end
label(wx - 1.0, 0, 'Air', 1, 0, '', 0, 3);
label(60, 0, 'Air', 1, 0, '', 0, 3);
label(xo, yo + (ro + rs)/2, 'Aluminum, 1100', 1, 0, 'Primary', 1, 4);
label(xo, yo, steel_name, 1, 0, '', 0, 4);
xm = wx - wall - od;
xn = ox + wall + od;
label( xm, 0, 'Coil_wire', 1, 0, 'Secondary',  N_limb, 5);
label(-xm, 0, 'Coil_wire', 1, 0, 'Secondary',  N_limb, 5);
label( xn, 0, 'Coil_wire', 1, 0, 'Secondary', -N_limb, 5);
label(-xn, 0, 'Coil_wire', 1, 0, 'Secondary', -N_limb, 5);

mi_makeABC(7, 120, 0, 0, 0);
end

function draw_poly(P)
  n = size(P, 1);
  for k = 1:n
    mi_addnode(P(k,1), P(k,2));
  end
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
