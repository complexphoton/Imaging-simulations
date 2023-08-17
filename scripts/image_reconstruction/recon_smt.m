function recon_smt(data_dir)

syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'z_f_air', 'dx', 'wavelength_list', 'FOV_before_windowing', ...
    'epsilon_in', 'epsilon_eff', 'W_image', 'L_image', 'dy_image', 'dz_image', 'noise_amp', 'n_jobs', 'NA');

% In the calculation of r, the profile at the reference plane has the same
% phase at y = (W-FOV_before_windowing)/2, rather than (W-W_image)/2. Here
% y = 0 is the interface between the scattering structure and PML. Thus,
% the image reconstructed along y_image = 0 corresponds to y = (W-FOV_before_windowing)/2.
% Since our interested region is (W-W_image)/2 ~ (W+W_image)/2, we need to shift
% the reconstruction range along y by (FOV_before_windowing-W_image)/2.
y_image = (dy_image/2:dy_image:W_image) + (FOV_before_windowing - W_image)/2;

epsilon_medium_recon = epsilon_eff;
z_image = dz_image/2:dz_image:L_image;
% Define the image reconstruction grid Z and Y.
[Z, Y] = meshgrid(z_image, y_image);

[ny_image, nz_image] = size(Z);
psi = 0;

% Reshape Z and Y to column vectors, which are required for the Flatiron
% nufft.
Y = Y(:); Z = Z(:);

n_wavelength = round(length(wavelength_list)/n_jobs);
fprintf('reconstructing the image: ');
for job_jd = 1:n_jobs
    % get the wavelength list computed from job_jd
    wavelength_sublist = wavelength_list(((job_jd-1)*n_wavelength+1):min(job_jd*n_wavelength, length(wavelength_list)));

    textprogressbar(job_jd, job_jd/n_jobs*100);
    % load the corresponding subset of hyperspectral r.
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', num2str(job_jd)), 'hyperspectral_R_angular');

    for i = 1:length(wavelength_sublist)

        k0dx = 2*pi/wavelength_sublist(i)*dx;

        % For SMT, epsilon_medium_recon = epsilon_eff. It is important to
        % use the effective index rather than medium index here. Otherwise,
        % the imaging depth will slightly decrease and the target locations
        % along z will shift.
        channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_in, epsilon_medium_recon);

        % Select channels within NA.
        idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
        kz = channels.L.kxdx_prop(idx_NA)/dx;
        n_prop_NA = length(kz);

        % Load R and add noise to mimic the experiment condition.
        R = hyperspectral_R_angular{i};
        R = R + noise_amp*sqrt(mean(abs(R).^2, 'all'))*(randn(size(R))+1j*randn(size(R)));

        % The R is calculated with reference plane @ z_f_air in air for
        % mimicing experimental conditions.
        % Here, we shift the reference plane of r from the reference plane to the
        % surface.
        % In principle, we should add the t of the sample interface to
        % correct for the sample-induced aberration. But in the
        % case of air-water interface with NA = 0.5, t is diagonal and almost
        % constant over angles. Thus the effect of not correcting this
        % refraction should be negligible.
        R = reshape(exp(1i*kz*z_f_air), [], 1).*R.*reshape(exp(1i*kz*z_f_air), 1, []);

        % Obtain wavevectors in effective index.
        % Becasue ky remains unchanged during the refraction at the air-background interface,
        % the wavevectors are selected such that ky_{bg}(i) = ky_{air}(i).
        idx_bg_NA = (1:n_prop_NA) + round((channels.R.N_prop-n_prop_NA)/2);
        ky_bg = channels.R.kydx_prop(idx_bg_NA)/dx;
        kz_bg = channels.R.kxdx_prop(idx_bg_NA)/dx;

        %% Compute \psi(\omega) = sum_{ba} R_{ba}(\omega) exp(i((ky_b-ky_a)y+(kz_b - kz_a)z)) by nufft.
        % SMT = |sum_{\omega} \psi(\omega)|^2;

        % fy_{ba} = ky_bg_NA(b)-ky_bg_NA(a); fz_{ba} = -kz_bg_NA(b) - kz_bg_NA(a).
        % Here, ky_bg_NA/kz_bg_NA is a row vector with length n_prop_NA so
        % fy and fz are n_prop_NA-by-n_prop_NA matrices by implicit
        % expansion.
        fy = ky_bg.' - ky_bg;
        fz = -kz_bg.' - kz_bg;

        % f = finufft2d3(x,y,c,isign,eps,s,t) computes
        % f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
        % Note x[j], y[j] and c[j] are column vectors. The returned f is also a column vector.
        Ahr = finufft2d3(single(fy(:)), single(fz(:)), single(R(:)), 1, 1e-2, single(Y), single(Z));

        % Reshape the column vector to the 2D image.
        Ahr = reshape(Ahr, ny_image, nz_image);

        psi = psi + Ahr;

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

% Save the intensity and phase separately since we only plot the image
% intensity.
I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image.
normalization_factor = mean(I, 'all');
I = I/normalization_factor;
save(fullfile(recon_img_dir, 'smt'), 'y_image', 'z_image', 'I', 'phase_profile');

end
