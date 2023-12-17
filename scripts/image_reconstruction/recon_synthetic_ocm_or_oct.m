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
load(syst_data_path, 'z_f_air', 'dx', 'FOV_before_windowing', ...
    'epsilon_in', 'epsilon_eff', 'W_image', 'L_image', 'dy_image', 'dz_image', ...
    'noise_amp', 'n_jobs', 'n_wavelengths_per_job', 'NA');

%% Specify the reconstruction grid
y_image = (dy_image/2:dy_image:W_image) - W_image/2;
z_image = (dz_image/2:dz_image:L_image);

%% Check if Flatiron nufft exists
use_finufft = exist('finufft2d3', 'file');
if use_finufft
    fprintf(['reconstructing the ', upper(recon_method),' image with Flatiron nufft: ']);
else
    fprintf(['reconstructing the ', upper(recon_method),' image with MATLAB nufft: ']);
end

if strcmpi(recon_method, 'ocm')
    NA_method = 0.5;
else
    NA_method = 0.04;
end

%% Reconstruct the image
psi = 0; % complex synthetic OCM/OCT image amplitude
for job_id = 1:n_jobs
    % Display a text progress bar.
    textprogressbar(job_id, job_id/n_jobs*100);

    % Load the hyperspectral angular R for one job.
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', num2str(job_id)), ...
        'hyperspectral_R_angular', 'wavelength_list');

    for wavelength_idx = 1:n_wavelengths_per_job
        k0dx = 2*pi/wavelength_list(wavelength_idx)*dx;

        channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_in, epsilon_eff);

        % Select wavevectors within NA.
        idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
        idx_r_NA = abs(channels.L.kydx_prop(idx_NA)/(k0dx*sqrt(epsilon_in))) <= NA_method;
        n_prop_NA = sum(idx_r_NA);

        % Add a complex Gaussian noise to R.
        R = hyperspectral_R_angular{wavelength_idx};
        R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*(randn(size(R))+1j*randn(size(R)));
        R = R(idx_r_NA, idx_r_NA);

        % Obtain wavevectors in effective index.
        % Becasue ky remains unchanged during the refraction at the air-background interface,
        % the wavevectors are selected such that ky_{bg}(i) = ky_{air}(i).
        idx_bg_NA = (1:n_prop_NA) + round((channels.R.N_prop-n_prop_NA)/2);
        ky_bg = channels.R.kydx_prop(idx_bg_NA)/dx;

        % Compute the synthetic confocal reflection coefficient
        % Rc(y, omega) = sum_{ba} R_{ba}(omega)*exp(i(ky_b-ky_a)y) by NUFFT.
        % I_{OCM/OCT}(y, z) = |sum_{omega} Rc(y, omega)*exp(-2i*omega*t_f(z))|^2,
        % where exp(-2i*omega*t_f(z)) is the time gating factor.

        % fy_{ba} = ky_bg(b) - ky_bg(a);
        % Here, ky_bg is a row vector with length n_prop so
        % fy is an n_prop-by-n_prop matrix by implicit expansion.
        fy = ky_bg.' - ky_bg;

        if use_finufft
            % Rc_omega = finufft1d3(fy, R, isign, eps, y) computes
            % Rc_omega(m) = sum_j R(j)*exp(isign*i*fy(j)*y(m)), for m = 1, ..., ny.
            % fy, R, y, and Rc_omega are vectors.
            % isign = +1 or -1 and eps is the tolerance.
            Rc_omega = finufft1d3(single(fy(:)), single(R(:)), 1, 1e-2, single(y_image(:)));
        else
            % Rc_omega = nufftn(R, -fy/(2pi), {y}) computes
            % Rc_omega(m) = sum_j R(j)*exp(i*fy(j)*y(m)), for m = 1, ..., ny.
            % R and fy are column vectors with the same length.
            % {y} is a cell array, where y is a row vector.
            Rc_omega = nufftn(single(R(:)), -single(fy(:)/(2*pi)), {single(y_image)});
        end

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

        kz_in = channels.L.kxdx_prop(round(end/2))/dx;
        kz_eff = channels.R.kxdx_prop(round(end/2))/dx;
        time_gating_factor = exp(-2i*(kz_eff*z_image-kz_in*z_f_air));

        psi = psi + time_gating_factor.*reshape(Rc_omega, [], 1);
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
save(fullfile(recon_img_dir, ['synthetic_', lower(recon_method), '.mat']), 'y_image', 'z_image', 'I', 'phase_profile');
end
