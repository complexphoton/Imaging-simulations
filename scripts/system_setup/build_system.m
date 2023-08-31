% This script performs the following jobs:
% - Generate the locations of targets and low-index-contrast scatterers; 
% - Generate the permittivity profile;
% - Generate the ordering sequence via METIS; 
% - Save all necessary data for later usage.

clear
% It is recommended to run this simulation on a cluster.
data_dir = fullfile('.', 'data', 'small_system');
data_dir = fullfile('.', 'data', 'large_system');
%data_dir = fullfile('.', 'data', 'large_system_no_weak');

sim_params_path = fullfile(data_dir, 'set_simulation_parameters.m');
if isfile(sim_params_path)
    run(sim_params_path)
else
    error('Cannot find set_simulation_parameters.m! Please define the simulation parameters first or check your data_dir.')
end

%% Define the parameter for Planck-taper window.
FOV_before_windowing = FOV/(1-p); % specific to Planck-taper window; for Tukey-window it becomes FOV/(1-0.5*r). 

% generate an equally-spaced-frequency wavelength list
freq_over_c_min = 1/wavelength_max;
freq_over_c_max = 1/wavelength_min;
freq_over_c_list = linspace(freq_over_c_min, freq_over_c_max, n_wavelengths);
wavelength_list = 1./freq_over_c_list;

% check if the source profile go into PML.
if FOV+2*z_f_bg*tan(asin(NA/sqrt(epsilon_medium))) > W 
    error(['The system width is too narrow and the source amplitude at largest angle' ...
        ' could be non-negligible inside PML.'])
end

%% Generate scatterer locations
fprintf('generating scatterer locations...\n');
scatterer_min_separation = [[min_sep_between_target, min_sep_between_target_and_weak]; ...
    [min_sep_between_target_and_weak, min_sep_between_weak]];
scatterer_loc_list = generate_scatterer_locs(W, L, dx, [target_density; weak_density], ...
    [target_diam; weak_diam], scatterer_min_separation, [random_seed_target; random_seed_weak]);
target_locs = scatterer_loc_list{1};
weak_locs = scatterer_loc_list{2};

%% Generate the permittivity profile
scatterer_locs = [target_locs; weak_locs];
N_target = length(target_locs); N_weak = length(weak_locs);
N_scatterers = N_target + N_weak;

scatterer_diam = zeros(N_scatterers, 1);
scatterer_diam(1:N_target) = target_diam; 
scatterer_diam(N_target+(1:N_weak)) = weak_diam;

epsilon_scatterer = zeros(N_scatterers, 1);
epsilon_scatterer(1:N_target) = epsilon_target;
epsilon_scatterer(N_target+(1:N_weak)) = epsilon_weak;

fprintf('generating the permittivity profile...\n');
epsilon = generate_permittivity_profile(W, L, dx, dx, epsilon_medium, scatterer_locs, scatterer_diam, epsilon_scatterer);

epsilon_eff = mean(epsilon, 'all');

%% Compute the depth scaling factor and focal depth in air
% Due to the air-background interface, a beam focused to z_f_bg in the background
% will focus to z_f_air in air. Here, we compute the depth scaling factor z_f_bg/z_f_air
% and z_f_air. The derivation of the depth scaling factor can be found in
% Supplementary Sec. II.
wavelength = 2/(1/wavelength_list(1)+1/wavelength_list(end));
k0dx = 2*pi/wavelength*dx;
channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_in, epsilon_eff);

% Select channels within NA.
idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
kz = channels.L.kxdx_prop(idx_NA)/dx;
idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
idx_NA_bg = abs(channels.R.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;

kzdx = channels.L.kxdx_prop(idx_NA);
kzdx_bg = channels.R.kxdx_prop(idx_NA_bg);
depth_scaling = mean(sin(kzdx_bg)./sin(kzdx)); 
z_f_air = z_f_bg/depth_scaling; % focal depth in air

%% Define the scanning range of the focal spots
y_image_start = (W + 2*PML.npixels*dx - W_image)/2;
y_image_end = (W + 2*PML.npixels*dx + W_image)/2;

y_image = y_image_start:dy_image:y_image_end;
y_image_ocm = y_image-PML.npixels*dx-W/2;

y_image = y_image_start:dy_image_oct:y_image_end;
y_image_oct = y_image-PML.npixels*dx-W/2;

%% Generate the METIS ordering and reuse it during the R computation. 
if produce_metis_ordering
    fprintf('computing orderings...\n');
    generate_metis_ordering(data_dir);
end

syst_data_path = fullfile(data_dir, 'system_data.mat');
fprintf('saving the system data ...\n');
save(syst_data_path, 'L', 'W', 'FOV', 'FOV_before_windowing', 'z_f_bg', 'z_f_air', ...
    'depth_scaling', 'dx','epsilon', 'epsilon_in', 'epsilon_medium', 'epsilon_eff', ...
    'epsilon_target', 'epsilon_weak', 'PML', 'wavelength_list', 'noise_amp', ...
    'NA', 'NA_oct', 'W_image', 'L_image', 'dy_image', 'dz_image', 'dz_image_rcm', ...
    'p', 'y_image_oct', 'y_image_ocm', ...
    'random_seed_target', 'random_seed_weak', 'target_locs', 'weak_locs', ...
    'target_density', 'weak_density', 'n_jobs', 'n_wavelengths_per_job')


