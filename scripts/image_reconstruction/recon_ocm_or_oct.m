function recon_ocm_or_oct(data_dir, recon_method)
%recon_ocm_or_oct Reconstruct optical coherence microscopy (OCM) or 
%optical coherence tomography (OCT).
%    recon_ocm_or_oct(data_dir, recon_method)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/spatial_R stores the
%    spatial reflection matrices.
%  recon_method (character array; required):
%    The reconstruction method. Available choices are (case-insensitive):
%      'OCM'
%      'OCT'
%
%  === Outputs ===
%  There is no output arguement. The following variables will be saved in
%  data_dir/reconstructed_images/ocm.mat or oct.mat.
%  I (numeric matrix, real, single-precision):
%    The OCM/OCT image intensity. The intensity is normalized such
%    that its average is one.
%  phase_profile (numeric matrix, real, single-precision):
%    The phase of the OCM/OCT image. The complex image can be obtained from
%    psi = sqrt(I).*exp(1j*phase_profile).
%  y_image, z_image (row vector, real):
%    Coordinates of the OCM/OCT image. 
%
% === Notes ===
% - The code reconstructs the image between y = [-W_image/2, W_image/2] and
% z = [0, L_image], where W_image and L_image is defined during the system
% setup.

if ~(strcmpi(recon_method, 'ocm')||strcmpi(recon_method, 'oct'))
    error("The reconstruction method must be either ocm or oct.")
end

%% Load the system data
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'L_image', 'z_f_air', 'dx', 'dz_image',  ...
    'epsilon_in', 'epsilon_eff', 'noise_amp', 'n_jobs', 'n_wavelengths_per_job');

z_image = dz_image/2:dz_image:L_image;

%% Reconstruct the image
fprintf(['reconstructing the ', upper(recon_method),' image: ']);
rng(0) % fix the random seed for the Gaussian noise
psi = 0; % complex OCM/OCT image amplitude
for job_id = 1:n_jobs
    % Display a text progress bar.
    textprogressbar(job_id, job_id/n_jobs*100);

    % Load the hyperspectral spatial R for one job.
    if strcmpi(recon_method, 'ocm')
        R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA');
    else
        R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'low_NA');
    end
    load(fullfile(R_dir, [num2str(job_id), '.mat']), 'hyperspectral_R_spatial_diag', ...
        'wavelength_list', 'y_image');

    for i = 1:n_wavelengths_per_job
        % Load the hyperspectral confocal reflection coefficients.
        Rc_omega = hyperspectral_R_spatial_diag{i};

        % Add a complex Gaussian noise to R.
        Rc_omega = Rc_omega + noise_amp*sqrt(mean(abs(Rc_omega).^2, 'all'))*(randn(size(Rc_omega))+1j*randn(size(Rc_omega)));

        % Obtain wavevectors at normal incidence in the air and effective 
        % index for computing the time gating factor.
        % It is important to use epsilon_eff, rather than epsilon_medium,
        % to obtain kz_bg. Otherwise, the target locations along z will be inaccurate/shifted.
        k0dx = 2*pi/wavelength_list(i)*dx;
        channels = mesti_build_channels(1, 'TM', 'periodic', k0dx, epsilon_in, epsilon_eff);
        kz_in = channels.L.kxdx_prop(round(end/2))/dx;
        kz_bg = channels.R.kxdx_prop(round(end/2))/dx;

        % Compute the time-gating factor.

        % The incident pulse is at the air-medium interface z = 0 at time t = -z^{air}_f/v^{air}_g,
        % where v^{air}_g is the group velocity in air, and it would focus in the absence of the medium to
        % a depth z^{air}_f at time t = 0. It takes time z_image/v_g for
        % the pulse to travel from z = 0 to z = z_image in the medium,
        % where v_g is the group velocity in the medium and z_image is the
        % imaging depth, thus the pulse arrives at z = z_image at 
        % t = t_f = z_image/v_g - z^{air}_f/v^{air}_g.

        % To account for the numerical dispersion of the group
        % velocity, we set the time-gating factor exp(-2i*omega*t_f)
        % to exp(-2i(k^{eff}_z*z_image - k^{air}_z*z^{air}_f)),
        % where k^{eff}_z and k^{air}_z are the longitudinal components of
        % the wave vectors at normal incidence in the medium and in the air respectively.
        time_gating_factor = single(exp(-2i*(kz_bg*z_image-kz_in*z_f_air)));

        psi = psi + time_gating_factor.*Rc_omega;
    end
end
fprintf('done\n');

% Obtain the image intensity and phase
I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image such that the averaged image intensity is 1
I = I/mean(I, 'all');

%% Save the image data
recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end
save(fullfile(recon_img_dir, [lower(recon_method), '.mat']), 'y_image', 'z_image', 'I', 'phase_profile');
end

