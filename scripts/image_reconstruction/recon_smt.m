function recon_smt(data_dir, correct_for_index_mismatch)
%recon_smt Reconstruct scattering matrix tomography (SMT) in reflection mode.
%    recon_smt(data_dir)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/angular_R stores the
%    angular reflection matrices.
%  correct_for_index_mismatch (logical scalar; optional, defaults to true):
%    Whether to correct for the index-mismatch aberration. 
%
%  === Outputs ===
%  There is no output arguement. The following variables will be saved in
%  data_dir/reconstructed_images/smt.mat.
%  I (numeric matrix, real, single-precision):
%    The SMT image intensity. The intensity is normalized such
%    that its average is one.
%  phase_profile (numeric matrix, real, single-precision):
%    The phase of the SMT image. The complex image can be obtained from
%    psi = sqrt(I).*exp(1j*phase_profile).
%  y_image, z_image (row vector, real):
%    Coordinates of the SMT image. 
%
% === Notes ===
% - The code reconstructs the image between y = [-W_image/2, W_image/2] and
% z = [0, L_image], where W_image and L_image is defined during the system
% setup.

arguments
    data_dir 
    correct_for_index_mismatch = true
end

%% Load the system data.
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'z_f_air', 'dx', 'FOV_before_windowing', ...
    'epsilon_eff', 'epsilon_in', 'W_image', 'L_image', 'dy_image', 'dz_image', ...
    'depth_scaling', 'noise_amp', 'n_jobs', 'n_wavelengths_per_job');

%% Specify the reconstruction grid.
y_image = (dy_image/2:dy_image:W_image) - W_image/2;
z_image = dz_image/2:dz_image:L_image;
ny = length(y_image); nz = length(z_image);

if correct_for_index_mismatch
    % To correct for the refractive index mismatch at the air-medium interface,
    % we use the effective index of the scattering structure to synthesize input/output beams.
    epsilon_recon = epsilon_eff;
    z_image_recon = z_image;
else
    % Use the permittivity in air to synthesize input/output
    % beams that are focused in air.
    epsilon_recon = epsilon_in;

    % A beam focused in air at z will focus at z*depth_scaling in the
    % medium due to the refraction at the air-medium interface. 
    % Without the index-mismatch correction, the synthesized beams focus in
    % air. To reconstruct the image at depth z_image, the beam should focus
    % at z_image/depth_scaling.
    z_image_recon = z_image/depth_scaling;
end
[Z, Y] = meshgrid(z_image_recon, y_image);

%% Determine which NUFFT function to use.
use_finufft = exist('finufft2d3', 'file'); % check if Flatiron NUFFT exists.
if use_finufft
    fprintf('reconstructing the SMT image with Flatiron NUFFT: ');
else
    fprintf('reconstructing the SMT image with MATLAB NUFFT: ');
end

%% Reconstruct the SMT image.
rng(0) % fix the random seed for the Gaussian noise
psi = 0; % complex SMT image amplitude
for job_id = 1:n_jobs
    % Display a text progress bar.
    textprogressbar(job_id, job_id/n_jobs*100); 

    % Load the hyperspectral angular R for one job.
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', [num2str(job_id), '.mat']), ...
        'hyperspectral_R_angular', 'wavelength_list', 'ky_list', 'kz_list');

    for i = 1:n_wavelengths_per_job
        % Obtain the reflection matrix R at one frequency.
        R = hyperspectral_R_angular{i}; % a complex square matrix
        
        % Obtain wavevectors of the inputs/outputs of R.
        ky_air = ky_list{i}; % a row vector
        kz_air = kz_list{i};

        % Add a complex Gaussian noise to R.
        R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*(randn(size(R))+1j*randn(size(R)));

        % Shift the reference plane of R from z = z_f_air in air to z = 0.
        % In principle, we should add the transmission matrix (T) of the sample 
        % interface to correct for the sample-induced aberration. But in the
        % case of air-water interface with NA = 0.5, T is diagonal and almost
        % constant over angles. Thus the effect of not correcting this
        % refraction should be negligible.
        R = reshape(exp(1i*kz_air*z_f_air), [], 1).*R.*reshape(exp(1i*kz_air*z_f_air), 1, []);

        % Obtain wavevectors (ky_recon, kz_recon) in epsilon_recon for synthesizing input/output
        % beams.
        % Note that ky remains unchanged during the refraction at the air-background
        % interface. Thus, ky_recon = ky_air for epsilon_recon =
        % epsilon_eff or epsilon_in. The corresponding kz_recon is
        % calculated from the finite-difference dispersion relation.
        % When epsilon_recon = epsilon_in, kz_recon is simply kz_air.
        k0dx = 2*pi/wavelength_list(i)*dx;
        channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_recon);
        n_prop = length(kz_air);
        idx = (1:n_prop) + round((channels.N_prop-n_prop)/2);
        ky_recon = ky_air; 
        kz_recon = channels.kxdx_prop(idx)/dx;

        % Compute psi(y, z, omega) = sum_{ba} (R_{ba}(omega) exp(i((ky_b-ky_a)y+(kz_b-kz_a)z))) by NUFFT.
        % I_{SMT}(y, z) = |sum_{omega} psi(y, z, omega)|^2

        % fy_{ba} = ky_recon(b) - ky_recon(a); fz_{ba} = -kz_recon(b) - kz_recon(a).
        % Here, ky_recon and kz_recon are row vectors with length n_prop so
        % fy and fz are n_prop-by-n_prop matrices by implicit
        % expansion.
        fy = ky_recon.' - ky_recon;
        fz = -kz_recon.' - kz_recon;

        if use_finufft
            % psi_omega = finufft2d3(fy, fz, R, isign, eps, Y, Z) computes
            % psi_omega(k) = sum_j R(j)*exp(isign*i(fy(j)*Y(k) + fz(j)*Z(k))), for k = 1, ..., ny*nz.
            % fy, fz, R, Y, Z, and psi_omega are column vectors.
            % isign = +1 or -1 and eps is the tolerance.
            psi_omega = finufft2d3(single(fy(:)), single(fz(:)), single(R(:)), 1, 1e-2, single(Y(:)), single(Z(:)));
        else
            % psi_omega = nufftn(R, -[fy/(2pi), fz/(2pi)], {y, z}) computes
            % psi_omega(m, n) = sum_j R(j)*exp(i(fy(j)*y(m) + fz(j)*z(n))), for m = 1, ..., ny and n = 1, ..., nz.
            % R, fy, and fz are column vectors with the same length.
            % {y, z} is a cell array, where y and z are row vectors.

            % Two caveats: 
            % 1. The query points {y, z} can be a matrix, where each column corresponds to Y(:) or Z(:),  
            % but the speed will be much slower as of R2023b.
            % 2. The psi_omega that nufftn() returns is a column vector.
            psi_omega = nufftn(single(R(:)), -single([fy(:)/(2*pi), fz(:)/(2*pi)]), {single(y_image), single(z_image_recon)});
        end
            % Reshape the column vector to a 2D image.
            psi_omega = reshape(psi_omega, ny, nz);

        if ~correct_for_index_mismatch
            % We need an additional time-gating factor for reconstructing SMT
            % without index-mismatch correction. Without the correction, 
            % the incident pulse would focus in the absence of the medium to 
            % a depth z_image_recon at time t = 0. Due to the refraction at the air-medium interface,
            % the pulse focuses to a deeper depth z_image at time t = t_f = z_image/v_g-z_image_recon/v^{air}_g, 
            % where v_g is the group velocity in the medium and v^{air}_g is the group velocity in air.
            % To align the temporal gate with the spatial gate, we include an additional time-gating factor 
            % exp(-2i*omega*t_f).

            % To account for the numerical dispersion of the group
            % velocity, we set the time-gating factor exp(-2i*omega*t_f)
            % to exp(-2i(k^{eff}_z*z_image - k^{air}_z*z_image_recon)),
            % where k^{eff}_z and k^{air}_z are are the longitudinal components of 
            % the wave vectors at normal incidence in the medium and in the air respectively. 
            channels_eff = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_eff);
            kz_eff_norm = channels_eff.kxdx_prop(round(end/2))/dx;
            kz_air_norm = kz_air(round(end/2));
            time_gating_factor = exp(-2i*(kz_eff_norm*z_image-kz_air_norm*z_image_recon));
            psi_omega = psi_omega.*time_gating_factor;
        end

        psi = psi + psi_omega;
    end
end
fprintf('done\n');

% Obtain the image intensity and phase.
I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image such that the averaged intensity is 1.
I = I/mean(I, 'all');

%% Save the image data.
recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end

if correct_for_index_mismatch
    save(fullfile(recon_img_dir, 'smt.mat'), 'y_image', 'z_image', 'I', 'phase_profile');
else
    save(fullfile(recon_img_dir, 'smt_no_correction.mat'), 'y_image', 'z_image', 'I', 'phase_profile');
end

end
