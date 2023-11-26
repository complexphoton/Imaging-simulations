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

% A beam focused in air at z will focus at z*depth_scaling in the
% medium due to the refraction at the air-medium interface. 
% Without the index-mismatch correction, the synthesized beams focus in
% air. To reconstruct the image at depth z_image, the beam should focus
% at z_image/depth_scaling.
z_image_recon = z_image/depth_scaling;
ny = length(y_image); nz = length(z_image);
[Z, Y] = meshgrid(z_image_recon, y_image);

%% Determine which NUFFT function to use.
use_finufft = exist('finufft2d3', 'file'); % check if Flatiron NUFFT exists.
if use_finufft
    fprintf('reconstructing the RCM image with Flatiron nufft...\n');
else
    fprintf('reconstructing the RCM image with MATLAB nufft...\n');
end

%% Reconstruct the RCM image.
% Load the R from the center frequency.
job_id = round(n_jobs/2);
load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', [num2str(job_id), '.mat']), ...
    'hyperspectral_R_angular', 'ky_list', 'kz_list');
ky = ky_list{1}; kz = kz_list{1};
R = hyperspectral_R_angular{1};

% Add a complex Gaussian noise to R.
R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*randn(size(R), 'like', 1j);

% Shift the reference plane of R from z = z_f_air in air to z = 0.
R = reshape(exp(1i*kz*z_f_air), [], 1).*R.*reshape(exp(1i*kz*z_f_air), 1, []);

% Compute \psi(\omega_c) = sum_{ba} R_{ba}(\omega_c)
% exp(i((ky_b-ky_a)y+(kz_b - kz_a)z)) by NUFFT.
% RCM = |\psi(\omega_c)|^2;

% fy_{ba} = ky(b) - ky(a); fz_{ba} = -kz(b) - kz(a).
% Here, ky and kz is a row vector with length n_prop so
% fy and fz are n_prop-by-n_prop matrices by implicit
% expansion.
% Here, we use ky and kz in air because in RCM the beam is focused in air 
% rather than water.
fy = ky.' - ky;
fz = -kz.' - kz;

if use_finufft
    % psi = finufft2d3(fy, fz, R, isign, eps, y, z) computes
    % psi(k) = sum_j R(j) exp(isign*i (y(k) fy(j) + z(k) fz(j))), for k = 1, ..., length(fy).
    % fy, fz, R, y, z and psi are column vectors.
    % isign = +1 or -1 and eps is the tolerance.
    psi = finufft2d3(single(fy(:)), single(fz(:)), single(R(:)), 1, 1e-2, single(Y(:)), single(Z(:)));
else
    % psi = nufftn(R, -[fy/(2pi), fz/(2pi)], {y, z}) computes
    % psi(m, n) = sum_j R(j) exp(i (y(m) fy(j) + z(n)
    % fz(j))), for m = 1, ..., ny and n = 1, ..., nz.
    % R, fy and fz are column vectors with the same length.
    % {y, z} is a cell array.

    % Two caveats:
    % 1. The query points {y, z} can be a matrix, where each column corresponds to Y(:) or Z(:),
    % but the speed will be much slower as of R2023b.
    % 2. The psi that nufftn() returns is a column vector.
    psi = nufftn(single(R), -single([fy(:)/(2*pi), fz(:)/(2*pi)]), {single(y_image), single(z_image_recon)});
end

% Reshape the column vector to the 2D image.
I = abs(reshape(psi, ny, nz)).^2;

% Normalize the image such that the averaged intensity is 1.
I = I/mean(I, 'all');

%% Save the image data.
recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end
save(fullfile(recon_img_dir, 'rcm.mat'), 'y_image', 'z_image', 'I');
end
