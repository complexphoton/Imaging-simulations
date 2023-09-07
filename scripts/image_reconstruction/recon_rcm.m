function recon_rcm(data_dir)
%recon_rcm Reconstruct reflectance confocal microscopy (RCM).
%    recon_rcm(data_dir)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data and data_dir/hyperspectral_reflection_matrices/angular_R stores the
%    angular reflection matrices.
%
%  === Outputs ===
%  There is no output arguement. The following variables will be saved in
%  data_dir/reconstructed_images/rcm.mat.
%  I (numeric matrix, real, single-precision):
%    The RCM image intensity. The intensity is normalized such
%    that its average is one.
%  y_image, z_image (row vector, real):
%    Coordinates of the RCM image. 
%
% === Notes ===
% - The code reconstructs the image between y = [-W_image/2, W_image/2] and
% z = [0, L_image], where W_image and L_image is defined during the system
% setup.

%% Load the system data
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'z_f_air', 'depth_scaling', 'W_image', 'L_image', ...
    'dy_image', 'dz_image_rcm', 'noise_amp', 'n_jobs');

%% Specify the reconstruction grid
y_image = (dy_image/2:dy_image:W_image) - W_image/2;
z_image = (dz_image_rcm/2:dz_image_rcm:L_image);
% The image needs to be scaled along the z-axis due to index mismatch
[Z, Y] = meshgrid(z_image/depth_scaling, y_image);

%% Reconstruct the image
% Load the R from the center frequency
job_id = round(n_jobs/2);
load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', [num2str(job_id), '.mat']), ...
    'hyperspectral_R_angular', 'ky_list', 'kz_list');
ky = ky_list{1}; kz = kz_list{1};
R = hyperspectral_R_angular{1};

% Add a complex Gaussian noise to R
R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*randn(size(R), 'like', 1j);

% Shift the reference plane of R from z = z_f_air in air to z = 0
R = reshape(exp(1i*kz*z_f_air), [], 1).*R.*reshape(exp(1i*kz*z_f_air), 1, []);

% Compute \psi(\omega_c) = sum_{ba} R_{ba}(\omega_c) exp(i((ky_b-ky_a)y+(kz_b - kz_a)z)) by nufft
% RCM = |\psi(\omega_c)|^2;

% fy_{ba} = ky(b) - ky(a); fz_{ba} = -kz(b) - kz(a).
% Here, ky and kz is a row vector with length n_prop so
% fy and fz are n_prop-by-n_prop matrices by implicit
% expansion.
% Here, we use ky and kz in air because in RCM the beam is focused in air 
% rather than water.
fy = ky.' - ky;
fz = -kz.' - kz;

% f = finufft2d3(x,y,c,isign,eps,s,t) computes
% f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
% Note x[j], y[j] and c[j] are column vectors. The returned f is also a column vector.
Ahr = finufft2d3(single(fy(:)), single(fz(:)), single(R(:)), 1, 1e-2, single(Y(:)), single(Z(:)));

% Reshape the column vector to the 2D image
I = abs(reshape(Ahr, size(Y))).^2;

% Normalize the image such that the averaged intensity is 1
I = I/mean(I, 'all');

%% Save the image data
recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end
save(fullfile(recon_img_dir, 'rcm.mat'), 'y_image', 'z_image', 'I');
end
