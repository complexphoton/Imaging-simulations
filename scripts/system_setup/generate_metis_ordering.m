function generate_metis_ordering(data_dir)
%generate_metis_ordering Generate METIS orderings for the angular and
%spatial R computation.
%    generate_metis_ordering(data_dir)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data.
%    
%  === Outputs ===
%  There is no output arguement. The following files will be saved
%  under data_dir/orderings
%  angular_R.mat (row vector, postive integer):
%    The ordering for angular R.
%  spatial_R_high_NA.mat (row vector, positive integer):
%    The ordering for high-NA spatial R.
%  spatial_R_low_NA.mat (row vector, positive integer):
%    The ordering for low-NA spatial R.
%
% === Notes ===

load(fullfile(data_dir, 'system_data.mat'), 'n_jobs')

% Generate the ordering at the last wavelength (corresponds to the maximum frequency).
% At lower frequency, the size of R will be smaller thus one needs
% to pad redundant inputs/outputs to reuse the ordering at the max
% frequency.
job_id = n_jobs;

% generate the ordering for angular R.
compute_angular_R(data_dir, job_id, true);

% generate the ordering for spatial R with high-NA.
compute_spatial_R(data_dir, job_id, true, true);

% generate the ordering for spatial R with low-NA.
compute_spatial_R(data_dir, job_id, false, true);
end


