% - Define all parameters for reflection matrix computation and image reconstruction;

%% Define parameters for building the scattering structure.
% For the small system, the system thickness L is 1/5 of the imaging
% range of the large system (300 µm).
W = 140; % system width (micron) 
L = 60; % system length (micron)
dx = 0.1; % discretization grid size (micron)
epsilon_in = 1.0; % permittivity on the incident side
epsilon_medium = 1.33^2; % permittivity of the background medium

% scatterer parameters (diameter, density, permittivity)
target_diam = 0.3; % diameter (micron)
weak_diam = 0.7;
target_density = 0.01; % number of scatterers per unit area (micron^2)
weak_density = 0.25;
epsilon_target = 2.5^2; % scatterer permittivity
epsilon_weak = 1.4^2;

% random seeds for generating scatterer locations
random_seed_target = 2;
random_seed_weak = 0;

% minimum separation between scatterers
min_sep_between_target = 10*target_diam; % min. separation between target scatterers
min_sep_between_weak = weak_diam; % min. separation between low-index-contrast scatterers
min_sep_between_target_and_weak = (weak_diam+target_diam)*0.5; % min. separation between target and low-index-contrast scatterers

%% Define parameters for the reflection matrix computation.
% wavelength range and number of wavelengths
wavelength_min = 0.7;
wavelength_max = 1.0;
% The number of wavelengths should be large enough such that the
% vertical artifacts from periodical pulses is beyond imaging
% range. We estimate n_wavelengths from the empircal formula 
% z0/n_eff > L_image, where z0 = n_wavelengths/(2*n_eff*(1/wavelength_min-1/wavelength_max)) is the
% pulse spacing at normal incident, n_eff = <epsilon> is the
% effective index and L_image is the imaging range.
n_wavelengths = 90;

FOV = 80; % field of view (micron)
% Planck-taper window width.
% Note the p here is equivalent to the epsilon that is conventionally used.
% We do not use epsilon here to avoid the confusion.
p = 0.1; 
z_f_bg = 30; % depth of the focal plane in the background (micron)

% numerical apertures
NA = 0.5;
NA_oct = 0.04; % numerical apertures for the OCT simulation

% PML parameters
PML.npixels = 40;
PML.power_sigma = 3;
PML.sigma_max_over_omega = 60;
PML.power_kappa = 4;
PML.kappa_max = 1.7000;
PML.power_alpha = 1;
PML.alpha_max_over_omega = 0.7000;

% On a cluster, each job computes the reflection matrix at several wavelengths. 
% On a local machine, there is only one job. You can set n_wavelengths_per_job = n_wavelengths.
n_wavelengths_per_job = n_wavelengths; % number of wavelengths per job
n_jobs = round(n_wavelengths/n_wavelengths_per_job); % number of jobs submitted to the cluster
produce_metis_ordering = false; % set it to true to produce METIS ordering (METIS required)

%% Define parameters for the image reconstruction.
% reconstruction region. We reconstruct the image for y = [-W_image/2, W_image/2]
% and z = [0, L_image].
W_image = 60; % unit: micron
L_image = 60; 

% reconstruction resolution
dy_image = 0.2; % reconstruction resolution along the y-axis (micron)
dz_image = 0.2; % reconstruction resolution along the z-axis (micron)
% OCM and RCM has bad transverse and axial resolutions, respectively. 
% We increase the corresponding resolutions to speed up 
% the computation and reduce the data size.
dy_image_oct = 0.4;
dz_image_rcm = 0.5;

% Add noise to R for mimicing the experiment condition.
% The noise amplitude is determined in experiments by measuring 
% the relative fluctuation of R of a white paper over time. 
% The real and imaginary parts of measured noise are described by a complex 
% normal distribution, whose standard deviation gives noise_amp.
noise_amp = 0.1;
