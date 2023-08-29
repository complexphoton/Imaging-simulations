function recon_synthetic_ocm_or_oct(data_dir, recon_method)
%recon_smt Reconstruct synthetic OCM or OCT from the angular reflection
%matrix
%    recon_synthetic_ocm_or_oct(data_dir, recon_method)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/angular_R stores the
%    angular reflection matrices.
%  recon_method (character array; required):
%    The reconstruction method. Available choices are (case-insensitive):
%      'OCM'
%      'OCT'
%
%  === Outputs ===
%  There is no output arguement. The following variables will be saved in
%  data_dir/reconstructed_images/synthetic_ocm.mat or synthetic_oct.mat.
%  I (numeric matrix, real, single-precision):
%    The synthetic OCM/OCT image intensity. The intensity is normalized such
%    that its average is one.
%  phase_profile (numeric matrix, real, single-precision):
%    The phase of the synthetic OCM/OCT image. The complex image can be obtained from
%    psi = sqrt(I).*exp(1j*phase_profile).
%  y_image, z_image (row vector, real):
%    Coordinates of the synthetic OCM/OCT image. 
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
load(syst_data_path, 'z_f_air', 'dx', 'wavelength_list', 'FOV_before_windowing', ...
    'epsilon_in', 'epsilon_eff', 'W_image', 'L_image', 'dy_image', 'dz_image', ...
    'noise_amp', 'n_jobs', 'n_wavelengths_per_job', 'NA');

%% Specify the reconstruction grid
% In the calculation of R, the input/output transverse profiles are
% exp(1j*ky*(y-FOV_before_windowing/2)), where y = 0 corresponds to the
% system center. The reconstruction region is located between y = [-W_image/2,
% W_image/2]. Thus we need to shift the y-coordinate by
% FOV_before_windowing/2.
y_image = (dy_image/2:dy_image:W_image) + (FOV_before_windowing - W_image)/2;
z_image = (dz_image/2:dz_image:L_image);

if strcmpi(recon_method, 'ocm')
    NA_method = 0.5;
else
    NA_method = 0.04;
end

%% Reconstruct the image
fprintf('reconstructing images: ');
psi = 0;
for job_id = 1:n_jobs
    textprogressbar(job_id, job_id/n_jobs*100);

    % Load the hyperspectral R.
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', num2str(job_id)), 'hyperspectral_R_angular');
    % Compute the wavelength list.
    wavelength_sublist = wavelength_list(((job_id-1)*n_wavelengths_per_job+1):min(job_id*n_wavelengths_per_job, length(wavelength_list)));

    for wavelength_idx = 1:length(wavelength_sublist)
        k0dx = 2*pi/wavelength_sublist(wavelength_idx)*dx;

        % It is important to use the effective index rather than medium index here. 
        % Otherwise, the imaging depth will slightly decrease and the target locations
        % along z will shift.
        channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_in, epsilon_eff);

        % select wavevectors within NA.
        idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
        idx_r_NA = abs(channels.L.kydx_prop(idx_NA)/(k0dx*sqrt(epsilon_in))) <= NA_method;
        n_prop_NA = sum(idx_r_NA);

        % add noise
        R = hyperspectral_R_angular{wavelength_idx};
        R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*(randn(size(R))+1j*randn(size(R)));
        R = R(idx_r_NA, idx_r_NA);

        % Obtain wavevectors in effective index.
        % Becasue ky remains unchanged during the refraction at the air-background interface,
        % the wavevectors are selected such that ky_{bg}(i) = ky_{air}(i).
        idx_bg_NA = (1:n_prop_NA) + round((channels.R.N_prop-n_prop_NA)/2);
        ky_bg = channels.R.kydx_prop(idx_bg_NA)/dx;

        fy = ky_bg.' - ky_bg;
        Ahr = finufft1d3(single(fy(:)), single(R(:)), 1, 1e-2, single(y_image));
        
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

        % Compute time-gated reflection fields.
        psi = psi + time_gating_factor.*reshape(Ahr, [], 1);
    end
end
fprintf('done\n');

recon_img_dir = fullfile(data_dir, 'reconstructed_images');

if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end

% Define the y-coordinate for plotting, which puts y = 0 along the system
% center.
y_image = (dy_image/2:dy_image:W_image) - W_image/2;

I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image such that the averaged image intensity is 1.
normalization_factor = mean(I, 'all');
I = I/normalization_factor;

save(fullfile(recon_img_dir, ['synthetic_', recon_method]), 'y_image', 'z_image', 'I', 'phase_profile');
end
