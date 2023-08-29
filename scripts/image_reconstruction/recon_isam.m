function recon_isam(data_dir)
%recon_smt Reconstruct interferometric synthetic aperture microscopy (ISAM).
%    recon_isam(data_dir)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/spatial_R stores the
%    angular reflection matrices.
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
load(syst_data_path, 'dx', 'wavelength_list', 'epsilon_in', 'epsilon_eff', ...
    'W_image', 'L_image', 'dy_image', 'dz_image', 'noise_amp', 'n_jobs', ...
    'n_wavelengths_per_job', 'z_f_air');

%% Specify the reconstruction grid
y_image = dy_image/2:dy_image:W_image;
z_image = dz_image/2:dz_image:L_image;
[Z, Y] = meshgrid(z_image, y_image);
[ny_image, nz_image] = size(Z);

fprintf('reconstructing images: ');
psi = 0;
for job_id = 1:n_jobs
    textprogressbar(job_id, job_id/n_jobs*100);

    % load the hyperspectral R
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA', num2str(job_id)), 'hyperspectral_R_spatial_diag');
    % compute the wavelength list
    wavelength_sublist = wavelength_list(((job_id-1)*n_wavelengths_per_job+1):min(job_id*n_wavelengths_per_job, length(wavelength_list)));

    for wavelength_idx = 1:length(wavelength_sublist)
        k0dx = 2*pi/wavelength_sublist(wavelength_idx)*dx;

        % add noise
        R_diag = hyperspectral_R_spatial_diag{wavelength_idx};
        R_diag = R_diag + noise_amp*sqrt(mean(abs(R_diag).^2, 'all'))*(randn(size(R_diag))+1j*randn(size(R_diag)));
        
        % transform the spatial R to angular R
        N = 2^(nextpow2(length(R_diag)));
        R_angular = fftshift(fft(R_diag, N));

        %% generate (Qy, Qz) in air
        Qy = (2*pi/N*(0:(N-1))-pi)/dy_image;

        % Compute kzdx_all for Qz from numerical dispersion
        sin_kzdx_over_two_sq = 0.25*(k0dx^2)*epsilon_in - sin((Qy/2)*dx/2).^2;
        kzdx_all = 2*asin(sqrt(sin_kzdx_over_two_sq));
        Qz = -2*kzdx_all/dx;
        
        % select propagating channels
        idx = imag(Qz)==0; 
        Qy = Qy(idx);
        Qz = Qz(idx);
        R_angular = R_angular(idx);
        
        R_angular = reshape(exp(-1i*Qz*z_f_air), [], 1).*R_angular;  

        %% generate (Qy, Qz) in the effective medium
        kydx_all = (Qy/2)*dx;
        k0dx2_epsilon = (k0dx^2)*epsilon_eff;
        sin_kzdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kydx_all/2).^2;
        kzdx_all = 2*asin(sqrt(sin_kzdx_over_two_sq));
        Qz_eff = -2*kzdx_all/dx;
        Qy_eff = Qy;
        %disp(length(Qy_eff(:))); disp(length(Qz_eff(:))); disp(length(R_angular(:)))
        % f = finufft2d3(x,y,c,isign,eps,s,t) computes 
        % f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
        % Note spatial frequencies x, y and fourier coefs c should be
        % column vectors. The returned f is also a column vector.
        Ahr = finufft2d3(single(Qy_eff(:)), single(Qz_eff(:)), single(R_angular(:)), 1, 1e-2, single(Y(:)), single(Z(:)));
        
        % Reshape the column vector to the 2D image.
        Ahr = reshape(Ahr, ny_image, nz_image);
        psi = psi + Ahr;
    end
end
fprintf('done\n');

y_image = (dy_image/2:dy_image:W_image) - W_image/2; % shift y = 0 to the system center.

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
save(fullfile(recon_img_dir, 'isam'), 'y_image', 'z_image', 'I', 'phase_profile');
end
