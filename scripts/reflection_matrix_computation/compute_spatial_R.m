function compute_spatial_R(data_dir, job_id, use_high_NA, produce_metis_ordering)
%compute_spatial_R Compute a subset of hyperspectral reflection matrix in spatial basis.
%    compute_spatial_R(data_dir, job_id, use_high_NA, produce_metis_ordering)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data.
%  job_id (positive integer; required; between 1 and total number of jobs):
%    Job ID, a number between 1 and n_jobs. It is used to determine which
%    wavelengths are computed for this job. On a local machine with one job, 
%    you can use job_id=1.
%  use_high_NA (logical scalar; required):
%    Whether to use high numerical aperture to compute R.
%  produce_metis_ordering (logical scalar; required):
%    A logical argument used for generating the METIS ordering. 
%    You should set it to false for the reflection matrix computation. 
%    Setting it to true will skip the factorization and only save the ordering data.
%    
%  === Outputs ===
%  There is no output arguement. The following variables will be saved
%  under data_dir/hyperspectral_reflection_matrices/spatial_R/[high_NA or low_NA]/[job_id].mat
%  hyperspectral_R_spatial_diag (cell array of column vectors, complex, single-precision):
%    The array of hyperspectral reflection matrices. Each element in the 
%    cell array corresponds to the diagonal of R at one wavelength.
%  wavelength_list (row vector, real):
%    Wavelengths for the R stored in hyperspectral_R_spatial_diag.
%  y_image (row vector, real):
%    Transverse scanning range of the focal spots. 
%
% === Notes ===

%% Load the system data.
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'epsilon', 'epsilon_medium', 'dx', 'epsilon_in', 'W', 'L', ...
    'z_f_air', 'PML', 'wavelength_list', 'n_jobs', 'n_wavelengths_per_job');

if job_id < 1 || job_id > n_jobs
    error('Invalid job id: job id should satisfy 1 <= job_id <= n_jobs')
end

fprintf('job id/total number of jobs = %d/%d.\n', job_id, n_jobs)

% load the focal spot scanning range for OCM/OCT
ordering_dir = fullfile(data_dir, 'orderings');
if use_high_NA
    load(syst_data_path, 'NA', 'y_image_ocm');
    R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'high_NA');
    ordering_path = fullfile(ordering_dir, 'spatial_R_high_NA.mat');
    y_image = y_image_ocm;
else
    load(syst_data_path, 'NA_oct', 'y_image_oct');
    R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'spatial_R', 'low_NA');
    ordering_path = fullfile(ordering_dir, 'spatial_R_low_NA.mat');
    y_image = y_image_oct;
    NA = NA_oct;
end

%% Specify options for producing METIS orderings (produce_metis_ordering = true)
%% or attempt to load orderings (produce_metis_ordering = false).
if produce_metis_ordering
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
    % Create the directory for saving the spatial R if not exist
    if ~exist(R_dir, 'dir')
        mkdir(R_dir);
    end
end

%% Specify the permittivity profile
nPML = PML.npixels; % number of PML pixels
ny_sca = round(W/dx); nz_sca = round(L/dx); % number of pixels of the scattering structure
ny_tot = ny_sca + 2*nPML;  % total number of pixels along the y-axis

% From left to right, we put nPML pixels of PML, one pixel of
% source/detection, nz_sca pixels of scattering medium and nPML pixels of
% PML. The total number of pixels along z axis is nz_sca + 2*nPML + 1;
% From top to bottom, we put nPML pixels of PML, ny_sca pixels of
% scattering medium and nPML pixels of PML. The total number of pixels
% along y axis is ny_sca + 2*nPML.
epsilon = [epsilon_in*ones(ny_tot, nPML+1),[epsilon_medium*ones(nPML, nz_sca);epsilon;epsilon_medium*ones(nPML, nz_sca)], epsilon_medium*ones(ny_tot, nPML)];

%% Specify the line source location
n_source = nPML + 1; % index of the source plane. The source is placed just in front of the PML.
z_source = -dx/2; % the source plane is half a pixel before the interface.

%% Generate input profiles on the focal plane
% Calculate Gaussian beam wasit for generating profile on the focal plane.
lambda_c = 2/(1/wavelength_list(1)+1/wavelength_list(end)); % center wavelength
w0 = lambda_c/(pi*NA); % Here, refractive index is already included in NA.

% Generate the list of E^in(z=z_f_air, y).
% Here, y.' is an ny-by-1 column vector, and y_f is a 1-by-M_in row vector.
% So, y.' - y_f is an ny-by-M_in matrix by implicit expansion.
% Then, E_f is an ny-by-M_in matrix whose m-th column is the cross section
% of the m-th Gaussian beam at z = z_f_air.
y = (0.5:ny_tot)*dx-(2*nPML*dx+W)/2;
E_f = exp(-(y.' - y_image).^2/(w0^2));

%% Specify the PML
syst.PML = PML;

%% Determine the wavelength range

syst.dx = dx;
opts.prefactor = -2i;

if produce_metis_ordering
    wavelength_list = wavelength_list(end);
    n_wavelengths = 1;

    fprintf('computing the ordering: ');
else
    wavelength_idx = (1:n_wavelengths_per_job)+(job_id-1)*n_wavelengths_per_job;
    wavelength_list = wavelength_list(wavelength_idx);
    n_wavelengths = length(wavelength_list); % number of wavelengths for this job

    fprintf('computing reflection matrices: ');
    hyperspectral_R_spatial_diag = cell(n_wavelengths, 1);
end

for i = 1:n_wavelengths
    textprogressbar(i, i/n_wavelengths*100);
    syst.wavelength = wavelength_list(i);

    %% Compute line source profiles
    k0dx = 2*pi/syst.wavelength*syst.dx;
    % Obtain wavevectors (ky, kz) in the incident medium within the NA.
    channels = mesti_build_channels(ny_tot, 'TM', 'periodic', k0dx, epsilon_in);
    idx_prop_NA = abs(channels.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
    kydx_NA = channels.kydx_prop(idx_prop_NA);
    kzdx_NA = channels.kxdx_prop(idx_prop_NA);
    
    % Obtain transverse profiles of propagating plane waves.
    % Each column of u is one transverse profile. Different columns are orthonormal.
    u = channels.fun_u(kydx_NA);

    % Project E^in(z_f_air, y) onto the propagating plane waves.
    E_f_prop = u'*E_f;

    % Backpropagate the profile from z = z_f_air to z = 0 to obtain the
    % surface source profile.
    E_s_prop = exp(1i*kzdx_NA(:)/syst.dx*(z_source-z_f_air)).*E_f_prop;

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

    %% Calculate D. 
    % For a homogeneous space, the length of the simulation domain doesn't 
    % matter, so we choose a minimal thickness of nz_temp = n_source +
    % nPML.
    syst.epsilon = epsilon_in*ones(ny_tot, n_source + nPML);
    opts.ordering = []; % clear the ordering for D.
    D = mesti(syst, B_struct, C, [], opts);

    %% Calculate R or the METIS ordering.
    if ~produce_metis_ordering
        opts.ordering = ordering; % specify the METIS ordering
    end
    syst.epsilon = epsilon;
    [R, info] = mesti(syst, B_struct, C, D, opts);

    if ~produce_metis_ordering
        % convert R to single-precision
        hyperspectral_R_spatial_diag{i} = single(diag(R));
    end
end
fprintf('done\n');

if produce_metis_ordering
    % save the ordering
    ordering = info.ordering;
    save(ordering_path, 'ordering');
else
    save(fullfile(R_dir, num2str(job_id)), 'hyperspectral_R_spatial_diag', ...
        'wavelength_list', 'y_image');
end

end

