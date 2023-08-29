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
load(syst_data_path, 'z_f_air', 'dx', 'wavelength_list', 'FOV_before_windowing', ...
    'depth_scaling', 'epsilon_in', 'W_image', 'L_image', 'dy_image', 'dz_image_rcm',...
    'noise_amp', 'n_jobs', 'n_wavelengths_per_job', 'NA');

%% Specify the reconstruction grid
% In the calculation of R, the input/output transverse profiles are
% exp(1j*ky*(y-FOV_before_windowing/2)), where y = 0 corresponds to the
% system center. The reconstruction region is located between y = [-W_image/2,
% W_image/2]. Thus we need to shift the y-coordinate by
% FOV_before_windowing/2.
y_image = (dy_image/2:dy_image:W_image) + (FOV_before_windowing - W_image)/2;
z_image = (dz_image_rcm/2:dz_image_rcm:L_image);
% The image needs to be scaled along the z-axis due to index mismatch
[Z, Y] = meshgrid(z_image/depth_scaling, y_image);
[ny_image, nz_image] = size(Z);

%% Reconstruct the image
% use the center frequency for reconstruction
job_id = round(n_jobs/2);
wavelength = wavelength_list((job_id-1)*n_wavelengths_per_job+1);
k0dx = 2*pi/wavelength*dx;

% load R and add noise.
load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', num2str(job_id)), 'hyperspectral_R_angular');
R = hyperspectral_R_angular{1};
R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*(randn(size(R))+1j*randn(size(R)));

% In RCM, the beam is focused in air rather than water for
% mimicing experimental conditions.
epsilon_medium_recon = epsilon_in;
channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_in, epsilon_medium_recon);

% select wavevectors within NA.
idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
kz = channels.L.kxdx_prop(idx_NA)/dx;
n_prop_NA = length(kz);

% The R is calculated with reference plane @ z_f_air in air for
% mimicing experimental conditions.
% Here, we shift the reference plane of R from the reference plane to the
% surface.
R = reshape(exp(1i*kz*z_f_air), [], 1).*R.*reshape(exp(1i*kz*z_f_air), 1, []);

% Obtain wavevectors in effective index.
% Becasue ky remains unchanged during the refraction at the air-background interface,
% the wavevectors are selected such that ky_{bg}(i) = ky_{air}(i).
idx_bg_NA = (1:n_prop_NA) + round((channels.R.N_prop-n_prop_NA)/2);
ky_bg = channels.R.kydx_prop(idx_bg_NA)/dx;
kz_bg = channels.R.kxdx_prop(idx_bg_NA)/dx;

% fy_{ba} = ky_bg_NA(b)-ky_bg_NA(a); fz_{ba} = -kz_bg_NA(b) - kz_bg_NA(a).
% Here, ky_bg_NA/kz_bg_NA is a row vector with length n_prop_NA so
% fy and fz are n_prop_NA-by-n_prop_NA matrices by implicit
% expansion.
fy = ky_bg.' - ky_bg;
fz = -kz_bg.' - kz_bg;

% f = finufft2d3(x,y,c,isign,eps,s,t) computes
% f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
% Note x[j], y[j] and c[j] are column vectors. The returned f is also a column vector.
Ahr = finufft2d3(single(fy(:)), single(fz(:)), single(R(:)), 1, 1e-2, single(Y(:)), single(Z(:)));

% reshape the column vector to the 2D image.
I = abs(reshape(Ahr, ny_image, nz_image)).^2;

recon_img_dir = fullfile(data_dir, 'reconstructed_images');
if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end

% Define the y-coordinate for plotting, which puts y = 0 along the system
% center.
y_image = (dy_image/2:dy_image:W_image) - W_image/2;

% Normalize the image such that the averaged intensity is 1.
normalization_factor = mean(I, 'all');
I = I/normalization_factor;
save(fullfile(recon_img_dir, 'rcm'), 'y_image', 'z_image', 'I');
end
