%% ------------------------------------------------------------------------
%  SETTINGS & DATA‐FOLDER LAYOUT
%  [1] data/
%       ├─ distances.mat        % contains mtx_euc_dis
%       ├─ structural.mat       % contains full_connectome_no_symm
%       └─ model_outputs/       % one .mat per subject, each with A_hat, A, TR, output.eff_conn.NoiseVar, etc.
%           ├─ sub001.mat
%           ├─ sub002.mat
%           └─ …
%  ------------------------------------------------------------------------

clearvars, close all; clc

dataDir = fullfile(pwd,'data');
outDir = fullfile(dataDir,'regressed_001_01_sim62131');

% get list of your subject‐model .mat files
d = dir(fullfile(outDir,'*.mat'));
files = {d.name};

% create an output folder for the A matrices
A_dir = fullfile(outDir, 'A_matrices');
if ~exist(A_dir, 'dir')
    mkdir(A_dir)
end

for iSub = 1:numel(files)
    % Load this subject's fitted A, TR, noise variance
    subj = load(fullfile(outDir,files{iSub}));
    A = subj.A;
    
    % Define a filename (e.g., sub001_A.csv)
    [~, name, ~] = fileparts(files{iSub});
    csvFile = fullfile(A_dir, [name '_A.csv']);
    
    % Save A as CSV
    writematrix(A, csvFile);
end

disp('All A matrices saved.');