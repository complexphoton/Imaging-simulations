function recon_isam(data_dir)
%recon_smt Reconstruct interferometric synthetic aperture microscopy (ISAM).
%    recon_isam(data_dir)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/spatial_R/
%    high_NA stores the high-NA spatial reflection matrices.
%
%  === Outputs ===
%  There is no output arguement. The following variables will be saved in
%  data_dir/reconstructed_images/isam.mat.
%  I (numeric matrix, real, single-precision):
%    The ISAM image intensity. The intensity is normalized such
%    that its average is one.
%  phase_profile (numeric matrix, real, single-precision):
%    The phase of the ISAM image. The complex image can be obtained from
%    psi = sqrt(I).*exp(1j*phase_profile).
%  y_image, z_image (row vector, real):
%    Coordinates of the ISAM image. 
%
% === Notes ===
% - The code reconstructs the image between y = [-W_image/2, W_image/2] and
% z = [0, L_image], where W_image and L_image is defined during the system
% setup.

%% Load the system data
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'dx', 'epsilon_in', 'epsilon_eff', ...
    'W_image', 'L_image', 'dy_image', 'dz_image', 'noise_amp', 'n_jobs', ...
    'n_wavelengths_per_job', 'z_f_air');

%% Specify the reconstruction grid
y_image = dy_image/2:dy_image:W_image;
z_image = dz_image/2:dz_image:L_image;
ny = length(y_image); nz = length(z_image);

%% Check if Flatiron nufft exists
use_finufft = exist('finufft2d3', 'file');
if use_finufft
    [Z, Y] = meshgrid(z_image, y_image);
    fprintf('reconstructing the ISAM image with Flatiron nufft: ');
else
    fprintf('reconstructing the ISAM image with MATLAB nufft: ');
end

rng(0) % fix the random seed for the Gaussian noise
psi = 0; % complex ISAM image amplitude
for job_id = 1:n_jobs
    % Display a text progress bar
    textprogressbar(job_id, job_id/n_jobs*100);

    % Load the hyperspectral spatial R for one job
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA', [num2str(job_id), '.mat']), ...
        'hyperspectral_R_spatial_diag', 'wavelength_list');

    for wavelength_idx = 1:n_wavelengths_per_job
        R_diag = hyperspectral_R_spatial_diag{wavelength_idx};

        % Add a complex Gaussian noise to R
        R_diag = R_diag + noise_amp*sqrt(mean(abs(R_diag).^2, 'all'))*(randn(size(R_diag))+1j*randn(size(R_diag)));
        
        % Transform R from the spatial domain to the angular domain by FFT
        N = 2^(nextpow2(length(R_diag)));
        R_angular = fftshift(fft(R_diag, N));

        %% Generate (Qy, Qz) in air
        Qy = (2*pi/N*(0:(N-1))-pi)/dy_image;

        % Compute Qz using numerical dispersion
        k0dx = 2*pi/wavelength_list(wavelength_idx)*dx;
        sin_kzdx_over_two_sq = 0.25*(k0dx^2)*epsilon_in - sin((Qy/2)*dx/2).^2;
        kzdx_all = 2*asin(sqrt(sin_kzdx_over_two_sq));
        Qz = -2*kzdx_all/dx;
        
        % Select propagating components
        idx = imag(Qz)==0; 
        Qy = Qy(idx); Qz = Qz(idx);
        R_angular = R_angular(idx);
        
        %% Shift the reference plane from z = z_f_air in air to z = 0
        R_angular = reshape(exp(-1i*Qz*z_f_air), [], 1).*R_angular;  

        %% Generate (Qy, Qz) in the effective medium
        Qy_eff = Qy;

        % Compute Qz_eff using numerical dispersion
        kydx_all = (Qy/2)*dx;
        k0dx2_epsilon = (k0dx^2)*epsilon_eff;
        sin_kzdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kydx_all/2).^2;
        kzdx_all = 2*asin(sqrt(sin_kzdx_over_two_sq));
        Qz_eff = -2*kzdx_all/dx;

        %% Compute \psi(\omega) = sum_{ba} R_{ba}(\omega) exp(i(Qy_{ba}y+Qz_{ba}z)) by nufft
        % ISAM = |sum_{\omega} \psi(\omega)|^2;

        if use_finufft
             % f = finufft2d3(x,y,c,isign,eps,s,t) computes
             % f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
             % Note x[j], y[j] and c[j] are column vectors. The returned f is also a column vector.
            psi_omega = finufft2d3(single(Qy_eff(:)), single(Qz_eff(:)), single(R_angular(:)), 1, 1e-2, single(Y(:)), single(Z(:)));
        else
            psi_omega = nufftn(single(R_angular), [-single(Qy_eff(:)/(2*pi)), -single(Qz_eff(:)/(2*pi))], {single(y_image), single(z_image)});
        end

        % Reshape the column vector to the 2D image
        psi_omega = reshape(psi_omega, ny, nz);
        psi = psi + psi_omega;
    end
end
fprintf('done\n');

% Obtain the image intensity and phase
I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image such that the averaged image intensity is 1
I = I/mean(I, 'all');

% Shift the center of y_image to y = 0
y_image = y_image - W_image/2;

%% Save the image data
recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end
save(fullfile(recon_img_dir, 'isam.mat'), 'y_image', 'z_image', 'I', 'phase_profile');
end
