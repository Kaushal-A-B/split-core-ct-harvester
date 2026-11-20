% Draws the nonlinear model from ct_params.m and saves it to models/ct_model_dc.fem, no solving.
% Run it again after changing ct_params.m.

ct_setup;
ct_params;
openfemm(1);
ct_draw_model(0, 0, 0, []);
mf = fullfile(modeldir, 'ct_model_dc.fem');
mi_saveas(mf);
mi_close;
closefemm;
fprintf('saved %s\n', mf);
