% Stage 2 model of the design report: linear core, one joint, 50 Hz harmonic solve.
% Kept for comparison with the final model. Run from the repo root after addpath('src').

clear; clc;
ct_setup;

c   = 11.66;   t  = 34.98;          % limb width, stack depth (mm)
hw  = 22.19;   bw = 26.68;          % window height, width (mm)
lg  = 0.43;                         % air gap (mm)
wall   = 2.0;
d_wire = 0.611; od = 0.631;
layers = 2;    flange = 1.0;
Ncoil  = 64;                        % turns per limb
I1rms  = 300;
dcond  = 14.15;

wx = bw/2;  wy = hw/2;
ox = wx + c; oy = wy + c;
build = layers*od;
cy = wy - flange;
xi1 = wx - wall - build;  xi2 = wx - wall;
xo1 = ox + wall;          xo2 = ox + wall + build;

openfemm(0); newdocument(0);
mi_probdef(50, 'millimeters', 'planar', 1e-8, t, 30);

mi_drawrectangle(-ox,-oy, ox, oy);
mi_drawrectangle(-wx,-wy, wx, wy);
mi_addnode(-lg/2, wy); mi_addnode(lg/2, wy);
mi_addnode(-lg/2, oy); mi_addnode(lg/2, oy);
mi_addsegment(-lg/2, wy, -lg/2, oy);
mi_addsegment( lg/2, wy,  lg/2, oy);

r = dcond/2;
mi_addnode(-r,0); mi_addnode(r,0);
mi_addarc(-r,0, r,0, 180, 2.5);
mi_addarc( r,0,-r,0, 180, 2.5);

mi_drawrectangle( xi1,-cy,  xi2, cy);
mi_drawrectangle(-xi2,-cy, -xi1, cy);
mi_drawrectangle( xo1,-cy,  xo2, cy);
mi_drawrectangle(-xo2,-cy, -xo1, cy);

% the gap is drawn explicitly, so the core takes the initial permeability
mi_addmaterial('CRGO', 29557, 29557, 0, 0, 0, 0.3, 0, 0.95, 0, 0, 0, 0, 0);
mi_addmaterial('Coil_wire', 1, 1, 0, 0, 58, 0, 0, 1, 0, 0, 0, 1, d_wire);
mi_getmaterial('Air');
mi_getmaterial('Aluminum, 1100');

% FEMM takes peak values for harmonic problems
mi_addcircprop('Primary',   I1rms*sqrt(2), 1);
mi_addcircprop('Secondary', 0, 1);

mi_addblocklabel((wx+ox)/2, 0);  mi_selectlabel((wx+ox)/2, 0);
mi_setblockprop('CRGO', 1, 0, '', 0, 1, 0);  mi_clearselected;

mi_addblocklabel(0, (wy+oy)/2);  mi_selectlabel(0, (wy+oy)/2);
mi_setblockprop('Air', 0, 0.05, '', 0, 2, 0); mi_clearselected;

mi_addblocklabel(0, 9.5);        mi_selectlabel(0, 9.5);
mi_setblockprop('Air', 1, 0, '', 0, 3, 0);    mi_clearselected;

mi_addblocklabel(60, 0);         mi_selectlabel(60, 0);
mi_setblockprop('Air', 1, 0, '', 0, 3, 0);    mi_clearselected;

mi_addblocklabel(0, 0);          mi_selectlabel(0, 0);
mi_setblockprop('Aluminum, 1100', 1, 0, 'Primary', 0, 4, 1); mi_clearselected;

xm = (xi1+xi2)/2; xn = (xo1+xo2)/2;
mi_addblocklabel( xm,0); mi_selectlabel( xm,0); mi_setblockprop('Coil_wire',1,0,'Secondary',0,5, Ncoil); mi_clearselected;
mi_addblocklabel(-xm,0); mi_selectlabel(-xm,0); mi_setblockprop('Coil_wire',1,0,'Secondary',0,5, Ncoil); mi_clearselected;
mi_addblocklabel( xn,0); mi_selectlabel( xn,0); mi_setblockprop('Coil_wire',1,0,'Secondary',0,5,-Ncoil); mi_clearselected;
mi_addblocklabel(-xn,0); mi_selectlabel(-xn,0); mi_setblockprop('Coil_wire',1,0,'Secondary',0,5,-Ncoil); mi_clearselected;

mi_makeABC(7, 120, 0, 0, 0);

mi_saveas(fullfile(modeldir, 'ct_model_v2.fem'));
mi_analyze(1); mi_loadsolution;

cp = mo_getcircuitproperties('Secondary');
lam2 = cp(3);
I1pk = I1rms*sqrt(2);  N2 = 2*Ncoil;  w = 2*pi*50;
Lm_femm = N2*abs(lam2)/I1pk;
Voc_femm = w*abs(lam2)/sqrt(2);

lm  = 2*(hw+bw+2*c)/1000;  Ac = c*t*1e-6;  mu0 = 4*pi*1e-7;
mue = 29557/(29557*(lg/1000/lm)+1);
Lm_des = N2^2*mu0*mue*Ac/lm;
Voc_des = (I1rms/N2)*w*Lm_des;

B = mo_getb((wx+ox)/2, 0);
fprintf('Lm  : FEMM %.2f mH   vs design %.2f mH\n', Lm_femm*1e3, Lm_des*1e3);
fprintf('Voc : FEMM %.2f V    vs design %.2f V (rms)\n', Voc_femm, Voc_des);
fprintf('|B| in limb (peak) = %.3f T\n', norm(abs(B)));

mo_close;
mi_close;
closefemm;
