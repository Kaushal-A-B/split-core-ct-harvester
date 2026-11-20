function run_study(skip_femm, use_saved)
% Runs the whole study. run_study(true) skips FEMM and reuses the .dat files in results/.
% run_study(false, true) opens the model saved by ct_draw_save for the sine runs and their check.

if nargin < 1, skip_femm = false; end
if nargin < 2, use_saved = false; end

if ~skip_femm
  step('ct_build_model');
  step('ct_sine_femm', 'I1rms', 300, 'use_saved_model', use_saved);
  step('ct_sine_femm', 'I1rms', 500, 'use_saved_model', use_saved);
end
step('ct_load_sim', 'I1rms', 300, 'Tw', 75);
step('ct_load_sim', 'I1rms', 500, 'Tw', 125);
step('ct_thermal');
step('ct_bmap_plot');
step('ct_analytic');
if ~skip_femm
  step('ct_sine_check', 'I1rms', 300, 'use_saved_model', use_saved);
  step('ct_sine_check', 'I1rms', 500, 'use_saved_model', use_saved);
  step('ct_leakage');
  step('ct_offset_study');
end
end

function step(name, varargin)
  for k = 1:2:numel(varargin)
    eval(sprintf('%s = %g;', varargin{k}, varargin{k+1}));
  end
  fprintf('--- %s\n', name);
  eval(name);
end
