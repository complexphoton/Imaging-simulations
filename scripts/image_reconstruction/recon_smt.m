function recon_smt(data_dir)
%recon_smt Reconstruct scattering matrix tomography (SMT) in reflection mode.
%    recon_smt(data_dir)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/angular_R stores the
%    angular reflection matrices.
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

%% Load the system data
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'z_f_air', 'dx', 'FOV_before_windowing', ...
    'epsilon_eff', 'W_image', 'L_image', 'dy_image', 'dz_image', ...
    'noise_amp', 'n_jobs', 'n_wavelengths_per_job');

%% Specify the reconstruction grid
y_image = (dy_image/2:dy_image:W_image) - W_image/2;
z_image = dz_image/2:dz_image:L_image;
ny = length(y_image); nz = length(z_image);

%% Check if Flatiron nufft exists
use_finufft = exist('finufft2d3', 'file');
if use_finufft
    [Z, Y] = meshgrid(z_image, y_image);
    fprintf('reconstructing the SMT image with Flatiron nufft: ');
else
    fprintf('reconstructing the SMT image with MATLAB nufft: ');
end

%% Reconstruct the image
psi = 0; % complex SMT image amplitude
for job_id = 1:n_jobs
    % Display a text progress bar
    textprogressbar(job_id, job_id/n_jobs*100); 

    % Load the hyperspectral angular R for one job
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', [num2str(job_id), '.mat']), ...
        'hyperspectral_R_angular', 'wavelength_list', 'kz_list');

    for i = 1:n_wavelengths_per_job
        R = hyperspectral_R_angular{i};

        % Add a complex Gaussian noise to R
        R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*randn(size(R), 'like', 1j);

        % Shift the reference plane of R from z = z_f_air in air to z = 0
        % In principle, we should add the transmission matrix (T) of the sample 
        % interface to correct for the sample-induced aberration. But in the
        % case of air-water interface with NA = 0.5, T is diagonal and almost
        % constant over angles. Thus the effect of not correcting this
        % refraction should be negligible.
        kz = kz_list{i};
        R = reshape(exp(1i*kz*z_f_air), [], 1).*R.*reshape(exp(1i*kz*z_f_air), 1, []);

        % Obtain wavevectors in the effective index 
        % It is important to use wavevectors in the effective index rather 
        % than medium index to reconstruct the SMT image. Otherwise, the imaging depth 
        % will slightly decrease and the target locations along z will shift.
        k0dx = 2*pi/wavelength_list(i)*dx;
        channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_eff);

        % Becasue ky remains unchanged during the refraction at the air-background interface,
        % the wavevectors are selected such that ky_{bg}(i) = ky_{air}(i).
        n_prop = length(kz);
        idx = (1:n_prop) + round((channels.N_prop-n_prop)/2);
        ky_bg = channels.kydx_prop(idx)/dx;
        kz_bg = channels.kxdx_prop(idx)/dx;

        % Compute \psi(\omega) = sum_{ba} R_{ba}(\omega) exp(i((ky_b-ky_a)y+(kz_b - kz_a)z)) by nufft
        % SMT = |sum_{\omega} \psi(\omega)|^2;

        % fy_{ba} = ky_bg(b) - ky_bg(a); fz_{ba} = -kz_bg(b) - kz_bg(a).
        % Here, ky_bg and kz_bg is a row vector with length n_prop so
        % fy and fz are n_prop-by-n_prop matrices by implicit
        % expansion.
        fy = ky_bg.' - ky_bg;
        fz = -kz_bg.' - kz_bg;

        if use_finufft
             % f = finufft2d3(x,y,c,isign,eps,s,t) computes
             % f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
             % Note x[j], y[j] and c[j] are column vectors. The returned f is also a column vector.
            Ahr = finufft2d3(single(fy(:)), single(fz(:)), single(R(:)), 1, 1e-2, single(Y(:)), single(Z(:)));
        else
            Ahr = nufftn(single(R), [-single(fy(:)/(2*pi)), -single(fz(:)/(2*pi))], {single(y_image), single(z_image)});
        end

        % Reshape the column vector to the 2D image
        Ahr = reshape(Ahr, ny, nz);
        psi = psi + Ahr;
    end
end
fprintf('done\n');

% Obtain the image intensity and phase
I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image such that the averaged intensity is 1
I = I/mean(I, 'all');

%% Save the image data
recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end
save(fullfile(recon_img_dir, 'smt.mat'), 'y_image', 'z_image', 'I', 'phase_profile');
end
