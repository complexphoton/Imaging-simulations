% - Define all parameters for reflection matrix computation and image reconstruction;

W = 600; % system width
L = 400; % system length
FOV = 400; % field of view
z_f_bg = 150; % depth of the focal plane in the background

%% Define the wavelength range
wavelength_min = 0.7;
wavelength_max = 1.0;
% The number of wavelengths should be large enough such that the
% vertical artifacts from periodical pulses is beyond imaging
% range. We use the empircal formula z0/n_eff > L_image, where z0 =
% n_wavelength/(2*n_eff*(1/wavelength_min-1/wavelength_max)) is the
% pulse spacing at normal incident, n_eff = <epsilon> is the
% effective index and L_image is the imaging range, to estimate
% n_wavelngth.
n_wavelength = 450;

%% Define numerical apertures
NA = 0.5;
NA_oct = 0.04;

%% Define permittivity
epsilon_in = 1.0;
epsilon_medium = 1.33^2;

%% Parameters for the spatial R
% FOV gets smaller at the depth away from the focal plane. For the
% calculation of spatial R and image reconstruction, one shall estimate
% the appropriate FOV based on the formalism in the small
% system below.
FOV_min = 300;

%% Define scatterer parameters (diameter, density, epsilon).
target_diam = 0.3; % unit: micron
weak_diam = 0.7;
target_density = 0.002; % unit: micron^2
weak_density = 0.25;
epsilon_target = 2.5^2;
epsilon_weak = 1.4^2;

%% Define minimum distance between scatterers.
min_dist_between_weak = weak_diam; 
min_dist_between_target = 10*target_diam;
min_dist_between_target_and_weak = (weak_diam+target_diam)*0.5;

%% Add noise to R for mimicing the experiment condition.
% the noise is determined in experiment by measuring the relative fluctuation of r
% of a white paper over time. The experimentally measured noise is well
% described by a complex Gaussian and this noise_amp is the std of
% the gaussian relative to the mean.
noise_amp = 0.1;

%% Define PML parameters.
PML.npixels = 40;
PML.power_sigma = 3;
PML.sigma_max_over_omega = 60;
PML.power_kappa = 4;
PML.kappa_max = 1.7000;
PML.power_alpha = 1;
PML.alpha_max_over_omega = 0.7000;

%% Define discretiztion grid size
dx = 0.1;

produce_metis_ordering = false;

%% Define image reconstruction region and resolution.
W_image = 300;
L_image = 300;
dy_image = 0.2;
dz_image = 0.2;
dy_image_oct = 2;
dz_image_rcm = 0.5;

n_jobs = round(n_wavelength/2); % number of jobs submitted to the cluster

%% Define the random seed for generating scatterer locations.
random_seed_target = 9;
random_seed_weak = 0;

p = 0.1; % Note the p here is equivalent to the epsilon used in the conventional Planck-taper window.
