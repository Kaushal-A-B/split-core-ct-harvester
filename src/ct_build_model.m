% Draws the CT, solves it at 50 Hz with a linear core and compares Lm and Voc with a hand calculation.

clear; clc;
ct_setup;
ct_params;

mu_lin = 42345;

openfemm(1);
ct_draw_model(0, 0, f, mu_lin);
mi_saveas(fullfile(modeldir, 'ct_model_ac.fem'));
mi_analyze(1);
mi_loadsolution;

cp = mo_getcircuitproperties('Secondary');
I1pk = I1_design*sqrt(2);
lam2 = abs(cp(3));
Lm_femm  = N2*lam2/I1pk;
Voc_femm = w*lam2/sqrt(2);
mo_close;
mi_close;
closefemm;

lm  = 2*(hw + bw + 2*c)*1e-3;
Ac  = c*t*1e-6*stack_f;
lg  = 2*gap*1e-3;
mue = mu_lin/(mu_lin*lg/lm + 1);
Lm_hand  = N2^2*mu0*mue*Ac/lm;
Voc_hand = (I1_design/N2)*w*Lm_hand;

fprintf('Winding: %d turns/layer, %d + %d turns per limb, R1 = %.3f ohm at 20 C, copper %.1f g\n', ...
        tpl, n_full, n_part, R1_20, m_cu*1e3);
fprintf('Lm  : FEMM %.2f mH   hand calc %.2f mH\n', Lm_femm*1e3, Lm_hand*1e3);
fprintf('Voc : FEMM %.2f V    hand calc %.2f V (rms) at %g A\n', Voc_femm, Voc_hand, I1_design);
