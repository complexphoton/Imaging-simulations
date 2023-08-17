function recon_synthetic_ocm_or_oct(recon_method, data_dir)

if ~(strcmpi(recon_method, 'ocm')||strcmpi(recon_method, 'oct'))
    error("The reconstruction method must be either ocm or oct.")
end

syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'z_f_air', 'dx', 'wavelength_list', 'FOV_before_windowing', ...
    'epsilon_in', 'epsilon_eff', 'W_image', 'L_image', 'dy_image', 'dz_image', 'noise_amp', 'n_jobs','NA');

% In the calculation of r, the profile at the reference plane has the same
% phase at y = (W-FOV_before_windowing)/2, rather than (W-W_image)/2. Here
% y = 0 is the interface between the scattering structure and PML. Thus,
% the image reconstructed along y_image = 0 corresponds to y = (W-FOV_before_windowing)/2.
% Since our interested region is (W-W_image)/2 ~ (W+W_image)/2, we need to shift
% the reconstruction range along y by (FOV_before_windowing-W_image)/2.
y_image = (dy_image/2:dy_image:W_image) + (FOV_before_windowing - W_image)/2;

if strcmpi(recon_method, 'ocm')
    NA_method = 0.5;
else
    NA_method = 0.04;
end

epsilon_medium_recon = epsilon_eff;
z_image = (dz_image/2:dz_image:L_image);

psi = 0;
n_wavelength = round(length(wavelength_list)/n_jobs);

fprintf('reconstructing images: ');
for job_id = 1:n_jobs
    % get the wavelength list computed from job_id
    wavelength_sublist = wavelength_list(((job_id-1)*n_wavelength+1):min(job_id*n_wavelength, length(wavelength_list)));
    textprogressbar(job_id, job_id/n_jobs*100);

    % load the corresponding subset of hyperspectral r.
    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R', num2str(job_id)), 'hyperspectral_R_angular');
    for wavelength_idx = 1:length(wavelength_sublist)

        k0dx = 2*pi/wavelength_sublist(wavelength_idx)*dx;

        % For SMT, epsilon_medium_recon = epsilon_eff. It is important to
        % use the effective index rather than medium index here. Otherwise,
        % the imaging depth will slightly decrease and the target locations
        % along z will shift.
        channels = mesti_build_channels(round(FOV_before_windowing/dx), 'TM', 'periodic', k0dx, epsilon_in, epsilon_medium_recon);

        % Select channels within NA.
        idx_NA = abs(channels.L.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
        idx_r_NA = abs(channels.L.kydx_prop(idx_NA)/(k0dx*sqrt(epsilon_in))) <= NA_method;
        n_prop_NA = sum(idx_r_NA);

        % Load r and add noise to mimic the experiment condition.
        r_angular = hyperspectral_R_angular{wavelength_idx};
        r_angular = r_angular + noise_amp*sqrt(mean(abs(r_angular).^2, 'all'))*(randn(size(r_angular))+1j*randn(size(r_angular)));
        r_angular = r_angular(idx_r_NA, idx_r_NA);

        % Obtain wavevectors in effective index.
        % Becasue ky remains unchanged during the refraction at the air-background interface,
        % the wavevectors are selected such that ky_{bg}(i) = ky_{air}(i).
        idx_bg_NA = (1:n_prop_NA) + round((channels.R.N_prop-n_prop_NA)/2);
        ky_bg = channels.R.kydx_prop(idx_bg_NA)/dx;

        fy = ky_bg.' - ky_bg;
        Ahr = finufft1d3(single(fy(:)), single(r_angular(:)), 1, 1e-2, single(y_image));

        kz_in = channels.L.kxdx_prop(round(end/2))/dx;
        kz_eff = channels.R.kxdx_prop(round(end/2))/dx;
        time_gating_factor = exp(-2i*(kz_eff*z_image-kz_in*z_f_air));

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
