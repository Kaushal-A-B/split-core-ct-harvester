function run_study(skip_femm)
% Runs the whole study in order. run_study(true) skips the FEMM steps and reuses the saved
% .dat files in results/.

if nargin < 1, skip_femm = false; end

if ~skip_femm
  step('ct_build_model');
  step('ct_sine_femm', 'I1rms', 300);
  step('ct_sine_femm', 'I1rms', 500);
end
step('ct_load_sim', 'I1rms', 300, 'Tw', 75);
step('ct_load_sim', 'I1rms', 500, 'Tw', 125);
step('ct_thermal');
step('ct_bmap_plot');
step('ct_analytic');
if ~skip_femm
  step('ct_sine_check', 'I1rms', 300);
  step('ct_sine_check', 'I1rms', 500);
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
