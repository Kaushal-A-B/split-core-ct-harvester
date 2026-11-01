% Design values and material data, shared by every script.

mu0 = 4*pi*1e-7;

% core geometry (mm)
c  = 12;          % limb width
t  = 34;          % stack depth
hw = 23;          % window height
bw = 27;          % window width
gap = 0.2;        % each of the two joints

% grain oriented M-4, 0.27 mm (ATI Allegheny Ludlum GOES data, 2012)
rho_fe  = 7650;       % kg/m3
stack_f = 0.965;
Bsat    = 2.035;      % T
lam_d   = 0.27;       % mm

% B-H curve: peak permeability from the datasheet table, interpolated from 1.6 T to 1.85 T,
% then a saturation approach above
B_mu  = (2:16)/10;
mu_pk = [25816 28860 31128 32887 34379 35799 37269 38816 40356 41661 42345 41839 39366 33921 24249];
H_lo = 0.02/(mu0*12050);
H_tab = B_mu./(mu0*mu_pk);
H16 = H_tab(end);  H800 = 796;
B_mid = [1.65 1.70 1.75 1.80 1.85];
H_mid = H16*(H800/H16).^((B_mid - 1.6)/0.25);
J800 = 1.850 - mu0*H800;
H_hi = [1000 1500 2500 4000 8000 16000 32000 64000 160000 320000];
B_hi = Bsat - (Bsat - J800)*(H800./H_hi) + mu0*H_hi;
bh = [0 0; 0.02 H_lo; B_mu(:) H_tab(:); B_mid(:) H_mid(:); B_hi(:) H_hi(:)];

% core loss, W/kg, from the 50 Hz and 60 Hz tables; P/f = a + b f
B_loss = 1.0:0.1:1.8;
P50 = [.338 .410 .488 .575 .675 .792 .935 1.144 1.386];
P60 = [.203 .246 .293 .345 .405 .475 .561 .686 .830]*2.20462;
b_e = (P60/60 - P50/50)/10;
a_h = P50/50 - 50*b_e;

k_fe  = 29.9;       % W/mK
cp_fe = 502;        % J/kgK

% secondary winding, AWG 22
N_limb = 64;  N2 = 2*N_limb;
d_bare = 0.644;     % mm
od     = 0.70;      % mm
wall   = 2.0;       % mm, bobbin wall
flange = 1.0;       % mm
rho_cu20 = 1.724e-8;        % ohm m
alpha_cu = 0.00393;         % 1/K
Lw = hw - 2*flange;
tpl = floor(Lw/od + 1e-9);
n_full = 2*tpl;
n_part = N_limb - n_full;
off = wall + ((1:3) - 0.5)*od;
mlt = 2*(c + t) + 2*pi*off;
len_limb = tpl*mlt(1) + tpl*mlt(2) + n_part*mlt(3);
A_cu = pi/4*d_bare^2;
R1_20 = rho_cu20*(2*len_limb*1e-3)/(A_cu*1e-6);
m_cu  = 8960*(2*len_limb*1e-3)*(A_cu*1e-6);

% line conductor, ACSR Dog
dcond  = 14.15;     % mm
dsteel = 4.71;      % mm

f = 50;  w = 2*pi*f;
I1_design = 300;  I1_high = 500;      % A rms
RL_report = 9.88;                      % ohm

% rectifier, charger, battery
Vf = 0.5;           % V, SS34 at about 3 A
rd = 0.03;          % ohm
eta_buck = 0.90;
Vbat = [6.0 7.4 8.4];        % V, 2S pack empty / nominal / full
Rbat = 0.12;        % ohm
Cbat = 3.0;         % Ah
