function recon_ocm_or_oct(recon_method, data_dir)
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
load(syst_data_path, 'L_image', 'z_f_air', 'dx', 'dz_image', 'wavelength_list', ...
    'epsilon_in', 'epsilon_eff', 'noise_amp', 'n_jobs', 'n_wavelengths_per_job');

%% Specify the coordinates
if strcmpi(recon_method, 'ocm')
    load(syst_data_path, 'y_image_ocm');
    y_image = y_image_ocm;
else
    load(syst_data_path, 'y_image_oct');
    y_image = y_image_oct;
end
z_image = dz_image/2:dz_image:L_image;

%% Reconstruct the image
fprintf('reconstructing images: ');
psi = 0;
for job_id = 1:n_jobs
    textprogressbar(job_id, job_id/n_jobs*100);

    % load the hyperspectral R
    if strcmpi(recon_method, 'ocm')
        load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA', num2str(job_id)), 'hyperspectral_R_spatial_diag');
    else
        load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'low_NA', num2str(job_id)), 'hyperspectral_R_spatial_diag');
    end
    % compute the wavelength list
    wavelength_sublist = wavelength_list(((job_id-1)*n_wavelengths_per_job+1):min(job_id*n_wavelengths_per_job, length(wavelength_list)));
 
    for i = 1:n_wavelengths_per_job
        k0dx = 2*pi/wavelength_sublist(i)*dx;

        % It is important to use epsilon_eff, rather than epsilon_medium,
        % to obtain kz_eff. Otherwise, the target locations along z will be inaccurate/shifted.
        channels = mesti_build_channels(1, 'TM', 'periodic', k0dx, epsilon_in, epsilon_eff);
        R_diag = hyperspectral_R_spatial_diag{i};

        % add noise 
        R_diag = R_diag + noise_amp*sqrt(mean(abs(R_diag).^2, 'all'))*(randn(size(R_diag))+1j*randn(size(R_diag)));

        % The time-gating factor = exp(-i*omega*t_gated), where t_gated = 2(z-z_f_bg)/v_g_mid
        % and v_g_mid is the group velocity at the center frequency in epsilon_medium.
        % In this expression, z-z_f_bg is the distance between the imaging
        % plane at depth z and the reference plane at depth z_f_bg. Light
        % travels a round trip (=2(z-z_f_bg)) before it is captured thus
        % the time-gated reflection light from depth z should be sampled at t =
        % 2(z-z_f_bg)/v_g_mid.

        % However, that expression does not account for the numerical dispersion in v_g. 
        % Thus, we replace omega/v_g_mid with kz at the normal incident, which is calculated
        % at the corresponding frequency.
        kz_in = channels.L.kxdx_prop(round(end/2))/dx;
        kz_eff = channels.R.kxdx_prop(round(end/2))/dx;
        time_gating_factor = exp(-2i*(kz_eff*z_image-kz_in*z_f_air));

        % compute time-gated reflection fields
        psi = psi + time_gating_factor.*R_diag;
    end
end
fprintf('done\n');

% Save the intensity and phase separately since we only plot the image intensity.
I = abs(psi).^2;
phase_profile = angle(psi);

recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end

% Normalize the image such that the averaged image intensity is 1.
normalization_factor = mean(I, 'all');
I = I/normalization_factor;
save(fullfile(recon_img_dir, recon_method), 'y_image', 'z_image', 'I', 'phase_profile');
end

