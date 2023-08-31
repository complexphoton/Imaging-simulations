% - Define all parameters for reflection matrix computation and image reconstruction;
% - Please refer to the large-system script for more comments.

%% Define parameters for building the scattering structure.
W = 600;
L = 400; 
dx = 0.1;
epsilon_in = 1.0;
epsilon_medium = 1.33^2;

% scatterer parameters (diameter, density, permittivity)
target_diam = 0.3;
weak_diam = 0.7;
target_density = 0.002;
% This system removes the low-index-contrast scatterers
weak_density = 0; 
epsilon_target = 2.5^2;
epsilon_weak = 1.4^2;

% random seeds for generating scatterer locations
random_seed_target = 9;
random_seed_weak = 0;

% minimum separation between scatterers
min_sep_between_target = 10*target_diam;
min_sep_between_weak = weak_diam;
min_sep_between_target_and_weak = (weak_diam+target_diam)*0.5;

%% Define parameters for the reflection matrix computation.
% wavelength range and number of wavelengths
wavelength_min = 0.7;
wavelength_max = 1.0;
n_wavelengths = 450;

% numerical apertures
NA = 0.5;
NA_oct = 0.04;

FOV = 400; % field of view (micron)
p = 0.1; % Planck-taper window width.
z_f_bg = 150; % depth of the focal plane in the background (micron)

% PML parameters
PML.npixels = 40;
PML.power_sigma = 3;
PML.sigma_max_over_omega = 60;
PML.power_kappa = 4;
PML.kappa_max = 1.7000;
PML.power_alpha = 1;
PML.alpha_max_over_omega = 0.7000;

n_wavelengths_per_job = 2;
n_jobs = round(n_wavelengths/n_wavelengths_per_job);
produce_metis_ordering = false;

%% Define parameters for the image reconstruction.
% reconstruction region
W_image = 300;
L_image = 300;

% reconstruction resolution
dy_image = 0.2;
dz_image = 0.2;
dy_image_oct = 2;
dz_image_rcm = 0.5;

% noise amplitude
noise_amp = 0.1;