function compute_angular_R(data_dir, job_id, produce_metis_ordering)
%compute_angular_R Compute hyperspectral reflection matrix in angular basis.
%    compute_angular_R(data_dir, job_id, produce_metis_ordering)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data.
%  job_id (positive integer; required; between 1 and total number of jobs):
%    Job ID, a number between 1 and n_jobs. It is used to determine which
%    wavelengths are computed for this job. On a local machine with one job, 
%    you can use job_id=1.
%  produce_metis_ordering (logical scalar; required):
%    A logical argument used for generating the METIS ordering. 
%    You should set it to false for the reflection matrix computation. 
%    Setting it to true will skip the factorization and only save the ordering data.
%    
%  === Outputs ===
%  There is no output arguement. The following variables will be saved
%  under data_dir/hyperspectral_reflection_matrices/angular_R/[job_id].mat
%  hyperspectral_R_angular (cell array of matrices, complex, single-precision):
%    The array of hyperspectral reflection matrices. Each element in the 
%    cell array corresponds to the R at one wavelength.
%  wavelength_list (row vector, real):
%    Wavelengths for the R stored in hyperspectral_R_angular.
%  ky_list, kz_list (cell array of row vectors, real):
%    Wavevectors for the R stored in hyperspectral_R_angular. 
%
% === Notes ===

%% Load the system data.
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'epsilon', 'epsilon_medium', 'dx', 'epsilon_in', 'W', 'L',...
    'z_f_air', 'p', 'NA', 'PML', 'wavelength_list', 'FOV_before_windowing', ...
    'n_jobs', 'n_wavelengths_per_job');

if job_id < 1 || job_id > n_jobs
    error('Invalid job id: job id should satisfy 1 <= job_id <= n_jobs')
end

fprintf('job id/total number of jobs = %d/%d.\n', job_id, n_jobs)

%% Specify options for producing METIS orderings (produce_metis_ordering = true)
%% or attempt to load orderings (produce_metis_ordering = false).
ordering_dir = fullfile(data_dir, 'orderings');
ordering_path = fullfile(ordering_dir, 'angular_R.mat');
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
    % Create the directory for saving the angular R if not exist
    R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R');
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

%% Specify the line source location and window profile
n_source = nPML+1; % index of the source plane. The source is placed just in front of the PML.
z_source = -dx/2; % the source plane is half a pixel before the interface.

% Generate the Planck-taper window that will be applied to the source profile.
% This is used to suppress the source amplitude inside the PML.
ny_fov = round(FOV_before_windowing/dx);
win_profile = planckwin(ny_fov, p);

%% Specify the PML
syst.PML = PML;

% Get the maximum number of inputs over all wavelengths; 
% we will need to pad additional channels for the other wavelengths to
% use the METIS ordering.
channels = mesti_build_channels(ny_fov, 'TM', 'periodic', sqrt(epsilon_in)*2*pi/wavelength_list(end)*dx, epsilon_in);
N_prop_max = round(NA*channels.N_prop);

syst.dx = dx;
opts.prefactor = -2i;

%% Determine the wavelength range
if produce_metis_ordering
    wavelength_list = wavelength_list(end);
    n_wavelengths = 1;

    fprintf('computing the ordering: ');
else
    wavelength_idx = (1:n_wavelengths_per_job)+(job_id-1)*n_wavelengths_per_job;
    wavelength_list = wavelength_list(wavelength_idx);
    n_wavelengths = length(wavelength_list); % number of wavelengths for this job

    fprintf('computing reflection matrices: ');
    hyperspectral_R_angular = cell(n_wavelengths, 1);
    ky_list = cell(n_wavelengths, 1);
    kz_list = cell(n_wavelengths, 1);
end

for i = 1:n_wavelengths
    textprogressbar(i, i/n_wavelengths*100);
    syst.wavelength = wavelength_list(i);

    %% Generate input profiles on the focal plane
    k0dx = 2*pi/syst.wavelength*syst.dx;
    channels_fov = mesti_build_channels(ny_fov, 'TM', 'periodic', k0dx, epsilon_in);

    % Obtain wavevectors (ky, kz) in the incident medium within the NA.
    idx_prop_NA = abs(channels_fov.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;   
    kydx_NA = channels_fov.kydx_prop(idx_prop_NA); 
    kzdx_NA = channels_fov.kxdx_prop(idx_prop_NA);
    N_in = length(kydx_NA);
    ky_list{i} = kydx_NA/dx; kz_list{i} = kzdx_NA/dx;

    % Obtain transverse profiles of propagating plane waves.
    % Each column of u is one transverse profile. Different columns are orthonormal.
    E_i_prof = channels_fov.fun_u(kydx_NA);

    % Generate the list of input profile on the focal plane:
    % E^in(z=z_f_air, y)
    E_i = zeros(ny_tot, N_prop_max);
    % Index range of the FOV on the whole system width W + 2*PML_thickness. 
    y_i = (dx/2:dx:FOV_before_windowing) - FOV_before_windowing/2 + (W+2*nPML*dx)/2;
    y_i_idx = round(y_i/dx);
    % shift the center of the transverse profile to y = 0
    E_i_prof = E_i_prof.*exp(-1j*kydx_NA/dx*FOV_before_windowing/2);
    E_i(y_i_idx, 1:N_in) = reshape(win_profile, [], 1).*E_i_prof;
    
    %% Compute line source profiles
    % Get properties of propagating channels in the free space with system
    % width = ny_tot
    channels = mesti_build_channels(ny_tot, 'TM', 'periodic', k0dx, epsilon_in);
    u = channels.fun_u(channels.kydx_prop);

    % Project E^in(z_f_air, y) onto the propagating channels.
    E_i_prop = u'*E_i;

    % Backpropagate the profile from z = z_f_air to z = 0 to obtain the
    % surface source profile.
    kz_prop = channels.kxdx_prop/syst.dx;
    E_s_prop = reshape(exp(1i*kz_prop*(z_source-z_f_air)), [], 1).*E_i_prop;
    
    % Determine the line sources.
    % In a closed geometry with no PML in y, a line source of
    % -2i*nu(a)*u(:,a) generates outgoing waves with transverse profile
    % u(:,a). With PML in y, this is not strictly true but is sufficiently
    % accurate since E^in(z=0,y) decays exponentially in y.
    nu = reshape(channels.sqrt_nu_prop, [], 1).^2;
    B_L = u*(nu.*E_s_prop);
    
    % Check the source inside PML is small;
    if any(mean([abs(B_L(1:nPML, :)); abs(B_L((1:nPML)+ny_sca, :))], 1) > 1e-2*mean(abs(B_L(nPML+(1:ny_sca), :)), 1))
        warning('The source amplitude inside PML might be non-negligible. Please check your B_L.');
    end

    % In mesti(), B_struct.pos = [m1, n1, h, w] specifies the position of a
    % block source, where (m1, n1) is the index of the smaller-(y,x) corner,
    % and (h, w) is the height and width of the block. Here, we put line
    % sources (w=1) at n1 = n_source that spans the whole width of the
    % simulation domain (m1=1, h=ny).
    B_struct.pos = [1, n_source, ny_tot, 1];

    % B_struct.data specifies the source profiles inside such block, with
    % B_struct.data(:, a) being the a-th source profile.
    B_struct.data = B_L;

    C = 'transpose(B)';

    %% Calculate D. 
    % For a homogeneous space, the length of the simulation domain doesn't 
    % matter, so we choose a minimal thickness of nz_temp = n_source +
    % nPML.
    syst.epsilon = epsilon_in*ones(ny_tot, n_source+nPML);
    opts.ordering = []; % clear the ordering for D.
    D = mesti(syst, B_struct, C, [], opts);

    %% Calculate R or the METIS ordering.
    if ~produce_metis_ordering
        opts.ordering = ordering; % specify the METIS ordering
    end
    syst.epsilon = epsilon;
    [R, info] = mesti(syst, B_struct, C, D, opts);

    if ~produce_metis_ordering
        R = R(1:N_in, 1:N_in); % remove the redundant channels;
        R = flipud(R); % the output order is reversed from C = 'transpose(B)'.

        % approximately exclude the nu in C from C = 'transpose(B)'.
        sqrt_nu_NA_fov = channels_fov.sqrt_nu_prop(idx_prop_NA);
        R = reshape(1./sqrt_nu_NA_fov.^2, [], 1).*R;

        % convert R to single-precision
        hyperspectral_R_angular{i} = single(R);
    end
end
fprintf('done\n');

if produce_metis_ordering
    % save the ordering
    ordering = info.ordering;
    save(ordering_path, 'ordering');
else
    save(fullfile(R_dir, num2str(job_id)), 'hyperspectral_R_angular', ...
        'wavelength_list', 'ky_list', 'kz_list'); 
end
end
