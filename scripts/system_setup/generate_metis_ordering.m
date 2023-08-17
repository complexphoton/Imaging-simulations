function generate_metis_ordering(data_dir)
% Generate the ordering at the last wavelength (corresponds to the maximum frequency).
% At lower frequency, the size of matrix K will be smaller thus one needs
% to pad redundant inputs/outputs to reuse the ordering at the max
% frequency.
load(fullfile(data_dir, 'system_data.mat'), 'wavelength_list')

% use the last wavelength in wavelength_list, which corresponds to
% the maximum frequency with the most inputs.
n_jobs = length(wavelength_list);
job_id = n_jobs; % run the last job, which computes the last wavelength.

% generate the ordering for angular r
compute_angular_R(job_id, n_jobs, data_dir, true);

% generate the ordering for spatial r with high-NA.
compute_spatial_R(job_id, n_jobs, data_dir, true, true);

% generate the ordering for spatial r with low-NA.
compute_spatial_R(job_id, n_jobs, data_dir, false, true);
end


