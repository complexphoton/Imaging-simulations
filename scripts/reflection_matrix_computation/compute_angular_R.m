function compute_angular_R(job_id, n_jobs, data_dir, analysis_only)
% Compute hyperspectral angular reflection matrix
% - n_jobs: the number of jobs to be submitted on the cluster. 
% Each job computes the reflection matrix at several frequencies. 
% For example, we set n_jobs=225 for the large system simulation of 
% computing two frequencies per job and 450 frequencies in total. 
% On a local machine, you can set n_jobs=1.
% - job_id: a number between 1 and n_jobs. It is used to obtain the specific 
% frequency range for this job. For example, in our large system simulation 
% setting job_id=1 computes reflection matrices at the first two frequencies. 
% On a local machine, you can set job_id=1.
% - data_dir: the directory where system_data.mat is saved. 
% The reflection matrix data will be saved under this directory as well.
% - analysis_only: a logical argument used for generating the METIS ordering. 
% You should set it to false for the reflection matrix computation. 
% Setting it to true will skip the factorization and only save the ordering data.

syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'epsilon', 'epsilon_medium', 'dx', 'epsilon_in', 'W', 'L', 'z_f_air', 'p', 'NA', 'PML', 'wavelength_list', 'FOV_before_windowing');

R_dir = fullfile(data_dir, 'hyperspectral_reflection_matrices', 'angular_R');
ordering_dir = fullfile(data_dir, 'orderings');
ordering_path = fullfile(ordering_dir, 'angular_R');

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

% Get the wavelength list for this job.
wavelength_max = wavelength_list(end);
n_wavelength = round(length(wavelength_list)/n_jobs);

% If n_wavelength is not dividable by n_tasks, we will assign the remaining
% wavelengths to the last task.
wavelength_idx = ((job_id-1)*n_wavelength+1):min(job_id*n_wavelength, length(wavelength_list));
wavelength_list = wavelength_list(wavelength_idx);
n_wavelength = length(wavelength_list); % number of wavelength for this job

fprintf('job id/total number of jobs = %d/%d.\n', job_id, n_jobs)

% Define the PML and the permittivity inside PML.
syst.PML = PML;
nPML = PML.npixels;
epsilon_top = epsilon_medium; 
epsilon_bottom = epsilon_medium;

ny_sca = round(W/dx); nz_sca = round(L/dx);
ny_tot = ny_sca + 2*nPML; 

% From left to right, we put nPML pixels of PML, one pixel of
% source/detection, nz_sca pixels of scattering medium and nPML pixels of
% PML. The total number of pixels along z axis is nz_sca + 2*nPML + 1;
% From top to bottom, we put nPML pixels of PML, ny_sca pixels of
% scattering medium and nPML pixels of PML. The total number of pixels
% along y axis is ny_sca + 2*nPML.
epsilon = [epsilon_in*ones(ny_tot, nPML+1),[epsilon_top*ones(nPML, nz_sca);epsilon;epsilon_bottom*ones(nPML, nz_sca)], epsilon_medium*ones(ny_tot, nPML)];

n_source = nPML+1; % index of the source plane. The source is placed just in front of the PML.
z_source = -dx/2; % the source plane is half a pixel before the interface.

hyperspectral_R_angular = cell(n_wavelength, 1);
syst.dx = dx;
opts.prefactor = -2i;

% Generate the Planck-taper window that will be applied to the profile.
% This is used to suppress the source amplitude inside the PML.
ny_fov = round(FOV_before_windowing/dx);
win_profile = planckwin(ny_fov, p);

% Get the maximum number of inputs over the wavelength range; 
% we will need to pad additional channels for the other wavelengths to
% reuse the ordering.
channels_max = mesti_build_channels(ny_fov, 'TM', 'periodic', sqrt(epsilon_in)*2*pi/wavelength_max*dx, epsilon_in);
N_prop_max = round(NA*channels_max.N_prop);

% Index range of the FOV on the whole system width W + 2*PML_thickness. 
% It will be used to assign the profile on the focal plane.
y_i = (dx/2:dx:FOV_before_windowing) + (W+2*nPML*dx-FOV_before_windowing)/2;
y_i_idx = round(y_i/dx);

if analysis_only
    fprintf('computing the ordering: ');
else
    fprintf('computing reflection matrices: ');
end

for i = 1:n_wavelength
    textprogressbar(i, i/n_wavelength*100);
    syst.wavelength = wavelength_list(i);

    % Get properties of propagating channels in the incident medium for
    % generating profile on the focal plane.
    k0dx = 2*pi/syst.wavelength*syst.dx;
    channels_fov = mesti_build_channels(ny_fov, 'TM', 'periodic', k0dx, epsilon_in);

    % Select the inputs within NA.
    idx_prop_NA_fov = abs(channels_fov.kydx_prop/(k0dx*sqrt(epsilon_in))) <= NA;
    
    % Transverse profiles of the propagating channels within NA. 
    % Each column of u is one transverse profile. Different columns are orthonormal.
    kydx_NA_fov =  channels_fov.kydx_prop(idx_prop_NA_fov);
    E_i_prof = channels_fov.fun_u(kydx_NA_fov); % used for generating the profile on the focal plane.
    N_in = length(kydx_NA_fov);

    % Generate the list of E^in(z=z_f_air, y).
    E_i = zeros(ny_tot, N_prop_max);
    E_i(y_i_idx, 1:N_in) = reshape(win_profile, [], 1).*E_i_prof;

    % Get properties of propagating channels in the free space for
    % projection and backpropagation.
    channels = mesti_build_channels(ny_tot, 'TM', 'periodic', k0dx, epsilon_in);
    u = channels.fun_u(channels.kydx_prop);

    % Project E^in(z_f_air, y) onto the propagating channels.
    E_i_prop = u'*E_i;

    % Backpropagate the profile from z = z_f_air to z = 0 to obtain the
    % surface source profile.
    kx = channels.kxdx_prop/syst.dx; % list of wave numbers
    E_s_prop = reshape(exp(1i*kx*(z_source-z_f_air)), [], 1).*E_i_prop;
    
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

    % Calculate D
    % For a homogeneous space, the length of the simulation domain doesn't
    % matter, so we choose a minimal thickness of nz_temp = n_source + nPML
    syst.epsilon = epsilon_in*ones(ny_tot, n_source+nPML);
    opts.ordering = []; % clear the ordering for D.
    D = mesti(syst, B_struct, C, [], opts);

    % Calculate r or the ordering.
    syst.epsilon = epsilon;
    % reuse ordering.
    if ~analysis_only
        opts.ordering = ordering;
    end
    [R, info] = mesti(syst, B_struct, C, D, opts);

    if ~analysis_only
        R = R(1:N_in, 1:N_in); % remove the redundant channels;
        R = flipud(R); % the output order is reversed from C = 'transpose(B)'.

        % Approximately exclude the nu in C from C = 'transpose(B)'.
        sqrt_nu_NA_fov = channels_fov.sqrt_nu_prop(idx_prop_NA_fov);
        R = reshape(1./sqrt_nu_NA_fov.^2, [], 1).*R;

        % save the single-precision data
        hyperspectral_R_angular{i} = single(R);
    end
end
fprintf('done\n');

if analysis_only
    % save the ordering
    ordering = info.ordering;
    save(ordering_path, 'ordering');
else
    save(fullfile(R_dir, num2str(job_id)), 'hyperspectral_R_angular'); % save the hyperspectral r.
end

end
