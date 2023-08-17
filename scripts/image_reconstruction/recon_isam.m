function recon_isam(data_dir)
% Reconstruct interferometric synthetic aperture microscopy

syst_data_path = fullfile(data_dir, 'system_data.mat');

load(syst_data_path, 'dx', 'wavelength_list', 'epsilon_in', 'epsilon_eff', ...
    'W_image', 'L_image', 'dy_image', 'dz_image', 'noise_amp', 'n_jobs', 'z_f_air');

% Define the image reconstruction region.
W_recon = W_image;
y_image_start = (W_recon - W_image)/2;
y_image_end = (W_recon + W_image)/2;
y_image = y_image_start:dy_image:y_image_end;
z_image = dz_image/2:dz_image:L_image;

% Define the image reconstruction grid Z and Y.
[Z, Y] = meshgrid(z_image, y_image);

[ny_image, nz_image] = size(Z);

dy_f = dy_image; 

psi = zeros(ny_image, nz_image);

% Reshape Z and Y to column vectors, which are required for the Flatiron
% nufft.
Y = Y(:); Z = Z(:);
n_wavelength = round(length(wavelength_list)/n_jobs);

fprintf('reconstructing images: ');
for job_id = 1:n_jobs
    % Get the wavelength list computed from job_id.
    wavelength_sublist = wavelength_list(((job_id-1)*n_wavelength+1):min(job_id*n_wavelength, length(wavelength_list)));
    textprogressbar(job_id, job_id/n_jobs*100);

    load(fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA', num2str(job_id)), 'hyperspectral_R_spatial_diag');

    for wavelength_idx = 1:length(wavelength_sublist)
        wavelength = wavelength_sublist(wavelength_idx);

        k0 = 2*pi/wavelength;
        k0dx = k0*dx;

        % Load r and add noise to mimic the experiment condition.
        R_diag = hyperspectral_R_spatial_diag{wavelength_idx};
        R_diag = R_diag + noise_amp*sqrt(mean(abs(R_diag).^2, 'all'))*(randn(size(R_diag))+1j*randn(size(R_diag)));
        
        %% fft impl.
        N = 2^(nextpow2(length(R_diag)));
        R_angular = fftshift(fft(R_diag, N));

        % generate the Qy list
        Qy = (2*pi/N*(0:(N-1))-pi)/dy_f;

        % compute kzdx_all for Qz
        kydx_all = (Qy/2)*dx;
        k0dx2_epsilon = (k0dx^2)*epsilon_in;
        sin_kzdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kydx_all/2).^2;
        kzdx_all = 2*asin(sqrt(sin_kzdx_over_two_sq));
        Qz = -2*kzdx_all/dx;
        
        % select propagating channels
        idx = imag(Qz)==0; 
        Qy = Qy(idx);
        Qz = Qz(idx);
        R_angular = R_angular(idx);
        
        R_angular = reshape(exp(-1i*Qz*z_f_air), [], 1).*R_angular;  
        
        fy = Qy; 

        % compute kzdx_all
        kydx_all = (Qy/2)*dx;
        k0dx2_epsilon = (k0dx^2)*epsilon_eff;
        sin_kzdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kydx_all/2).^2;
        kzdx_all = 2*asin(sqrt(sin_kzdx_over_two_sq));
        fz = -2*kzdx_all/dx;
    
        % f = finufft2d3(x,y,c,isign,eps,s,t) computes 
        % f[k] = sum_j c[j] exp(+-i (s[k] x[j] + t[k] y[j])), for k = 1, ..., nk.
        % Note spatial frequencies x, y and fourier coefs c should be
        % column vectors. The returned f is also a column vector.
        Ahr = finufft2d3(single(fy(:)), single(fz(:)), single(R_angular(:)), 1, 1e-2, single(Y), single(Z));
        
        % Reshape the column vector to the 2D image.
        Ahr = reshape(Ahr, ny_image, nz_image);

        psi = psi + Ahr;
    end
end
fprintf('done\n');

y_image = (dy_image/2:dy_image:W_image) - W_image/2; % shift y = 0 to the system center.

I = abs(psi).^2;
phase_profile = angle(psi);

% Normalize the image such that the averaged image intensity is 1.
normalization_factor = mean(I, 'all');
I = I/normalization_factor;

recon_img_dir = fullfile(data_dir, 'reconstructed_images');

if ~exist(recon_img_dir, 'dir')
    mkdir(recon_img_dir);
end

save(fullfile(recon_img_dir, 'isam'), 'y_image', 'z_image', 'I', 'phase_profile');
end
