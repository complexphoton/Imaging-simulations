function compute_spatial_R(job_id, n_jobs, data_dir, use_high_NA, analysis_only)
%Compute hyperspectral spatial reflection matrix
% - job_id = 1, 2, ..., n_jobs: used for

% running multiple jobs in parallel on the cluster.
% Each job will compute a subset of the wavelength list based on its job id.
% - n_jobs: total number of jobs. It should be an integer between 1 and
% number of wavelengths.
% - data_dir: root of the data directory.
% - analysis_only: set it to true in the initialization script to obtain
% ordering without doing the factorization. set it to false for the r
% computation

syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'epsilon', 'epsilon_medium', 'dx', 'epsilon_in', 'z_f_air', 'PML', 'wavelength_list');

ordering_dir = fullfile(data_dir, 'orderings');
if use_high_NA
    load(syst_data_path, 'NA', 'y_image_ocm');
    R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA');
    ordering_path = fullfile(ordering_dir, 'spatial_R_high_NA');
    y_image = y_image_ocm;
else
    load(syst_data_path, 'NA_oct', 'y_image_oct');
    R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'low_NA');
    ordering_path = fullfile(ordering_dir, 'spatial_R_low_NA');
    y_image = y_image_oct;
    NA = NA_oct;
end

if analysis_only
    % Create the directory for saving the ordering if not exist
    if ~exist(ordering_dir, 'dir')
        mkdir(ordering_dir);
    end
    opts.store_ordering = true;
    opts.analysis_only = true;
    opts.use_METIS = true;
else
    if isfile(ordering_path)
        load(ordering_path, 'ordering')
    else
        warning('The ordering file is not found! Will use the default AMD ordering.')
        ordering = [];
    end
    % Create the directory for saving the angular R if not exist
    if ~exist(R_dir, 'dir')
        mkdir(R_dir);
    end
end

if job_id < 1 || job_id > n_jobs
    error('Invalid job id: job id should satisfy 1 <= job_id <= n_jobs')
end

% get the wavelength list for this job.
n_wavelength = round(length(wavelength_list)/n_jobs); % number of wavelength for this job
% if n_wavelength is not dividable by n_jobs, we will assign the remaining
% wavelengths to the last job.
wavelength_idx = ((job_id-1)*n_wavelength+1):min(job_id*n_wavelength, length(wavelength_list));

wavelength_sublist = wavelength_list(wavelength_idx);
n_wavelength = length(wavelength_sublist); % number of wavelength for this job

fprintf('job id/total number of jobs = %d/%d.\n', job_id, n_jobs)

% define the PML and the permittivity inside PML.
syst.PML = PML;
nPML = PML.npixels;
epsilon_top = epsilon_medium;
epsilon_bottom = epsilon_medium;

syst.dx = dx;
[ny_sca, nz_sca] = size(epsilon);
ny_tot = ny_sca + 2*nPML;

% From left to right, we put nPML pixels of PML, one pixel of
% source/detection, nz_sca pixels of scattering medium and nPML pixels of
% PML. The total number of pixels along z axis is nz_sca + 2*nPML + 1;
% From top to bottom, we put nPML pixels of PML, ny_sca pixels of
% scattering medium and nPML pixels of PML. The total number of pixels
% along y axis is ny_sca + 2*nPML.
epsilon = [epsilon_in*ones(ny_tot, nPML+1),[epsilon_top*ones(nPML, nz_sca);epsilon;epsilon_bottom*ones(nPML, nz_sca)], epsilon_medium*ones(ny_tot, nPML)];

n_source = nPML + 1; % index of the source plane. The source is placed just in front of the PML.
z_source = -dx/2; % the source plane is half a pixel before the interface.

lambda_c = 2/(1/wavelength_list(1)+1/wavelength_list(end)); % center wavelength

y = (0.5:ny_tot)*syst.dx;

% Calculate Gaussian beam wasit for generating profile on the focal plane.
w0 = lambda_c/(pi*NA); % Here, refractive index is already included in NA.

% Generate the list of E^in(z=z_f_air, y).
% Here, y.' is an ny-by-1 column vector, and y_f is a 1-by-M_in row vector.
% So, y.' - y_f is an ny-by-M_in matrix by implicit expansion.
% Then, E_f is an ny-by-M_in matrix whose m-th column is the cross section
% of the m-th Gaussian beam at z = z_f_air.
E_f = exp(-(y.' - y_image).^2/(w0^2));

hyperspectral_R_spatial_diag = cell(n_wavelength, 1);

opts.prefactor = -2i;

if analysis_only
    fprintf('computing the ordering: ');
else
    fprintf('computing reflection matrices: ');
end

for i = 1:n_wavelength
    textprogressbar(i, i/n_wavelength*100);
    syst.wavelength = wavelength_sublist(i);

    k0dx = 2*pi/syst.wavelength*syst.dx;
    % Get properties of propagating channels in the free space.
    channels = mesti_build_channels(ny_tot, 'TM', 'periodic', k0dx, epsilon_in);
    
    % Select the inputs within NA.
    idx_prop_NA = abs(channels.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
    
    % Transverse profiles of the propagating channels. Each column of u is
    % one transverse profile. Different columns are orthonormal.
    u = channels.fun_u(channels.kydx_prop(idx_prop_NA));

    % Project E^in(z_f_air, y) onto the propagating channels.
    E_f_prop = u'*E_f;

    % Backpropagate the profile from z = z_f_air to z = 0 to obtain the
    % surface source profile.
    kx = reshape(channels.kxdx_prop(idx_prop_NA)/syst.dx, [], 1); % list of wave numbers
    E_s_prop = exp(1i*kx*(z_source-z_f_air)).*E_f_prop;

    % Determine the line sources.
    % In a closed geometry with no PML in y, a line source of
    % -2i*nu(a)*u(:,a) generates outgoing waves with transverse profile
    % u(:,a). With PML in y, this is not strictly true but is sufficiently
    % accurate since E^in(z=0,y) decays exponentially in y.
    nu = reshape(channels.sqrt_nu_prop(idx_prop_NA), [], 1).^2;
    B_L = u*(nu.*E_s_prop);

    % In mesti(), B_struct.pos = [m1, n1, h, w] specifies the position of a
    % block source, where (m1, n1) is the index of the smaller-(y,x) corner,
    % and (h, w) is the height and width of the block. Here, we put line
    % sources (w=1) at n1 = n_source that spans the whole width of the
    % simulation domain (m1=1, h=ny).
    B_struct.pos = [1, n_source, ny_tot, 1];

    % B_struct.data specifies the source profiles inside such block, with
    % B_struct.data(:, a) being the a-th source profile.
    B_struct.data = B_L;

    % Check the source inside PML is small;
    if any(mean([abs(B_L(1:nPML, :)); abs(B_L((1:nPML)+ny_sca, :))], 1) > 1e-2*mean(abs(B_L(nPML+(1:ny_sca), :)), 1))
        warning('The source amplitude inside PML might be non-negligible. Please check your B_L.');
    end

    C = 'transpose(B)';

    % Calculate D. For a homogeneous space, the length of the simulation domain doesn't
    % matter, so we choose a minimal thickness of nz_temp = n_source + nPML
    syst.epsilon = epsilon_in*ones(ny_tot, n_source + nPML);
    opts.ordering = []; % clear the ordering for D.
    D = mesti(syst, B_struct, C, [], opts);

    % Calculate r or the ordering
    syst.epsilon = epsilon;

    % reuse ordering if not analysis_only mode
    if ~analysis_only
        opts.ordering = ordering;
    end

    [R, info] = mesti(syst, B_struct, C, D, opts);

    if ~analysis_only
        % save the single-precision data
        hyperspectral_R_spatial_diag{i} = single(diag(R));
    end
end
fprintf('done\n');

if analysis_only
    % save the ordering
    ordering = info.ordering;
    save(ordering_path, 'ordering');
else
    save(fullfile(R_dir, num2str(job_id)), 'hyperspectral_R_spatial_diag'); % save the hyperspectral r.
end

end

