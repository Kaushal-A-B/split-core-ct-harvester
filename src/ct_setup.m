% Paths shared by all scripts. Outputs go to results/, FEMM models to models/.
src = fileparts(mfilename('fullpath'));
root = fileparts(src);
outdir = fullfile(root, 'results');
modeldir = fullfile(root, 'models');

femm_path = getenv('FEMM_MFILES');
if isempty(femm_path), femm_path = 'C:\femm42\mfiles'; end
if exist(femm_path, 'dir'), addpath(femm_path); end

addpath(src);
if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(modeldir, 'dir'), mkdir(modeldir); end
cd(outdir);
