function compute_field_profile(data_dir, incident_profile_type)
%compute_field_profile Compute the field profile with the angular or
%spatial input.
%    compute_field_profile(data_dir, incident_profile_type)
%
%  === Input Arguments ===
%  data_dir (character array; required): 
%    The data directory where data_dir/system_data.mat stores
%    system data.
%  incident_profile_type (character array; required):
%    The incident profile type. Available choices are (case-insensitive):
%        'angular'
%        'spatial'
%    
%  === Outputs ===
%  There is no output arguement. The following variables will be saved
%  under data_dir/field_profile_[incident_profile_type].mat
%  I (numeric matrix, real, single-precision):
%    The field profile intensity. The intensity is normalized such
%    that its peak is one.
%  L_air (scalar, real):
%    The thickness of air on the left. The unit is micron.
%
% === Notes ===

if ~(strcmpi(incident_profile_type, 'angular')||strcmpi(incident_profile_type, 'spatial'))
    error("The incident profile type must be either angular or spatial")
end

%% Load the system data
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'W', 'L', 'epsilon', 'dx', 'PML', 'NA', ...
    'z_f_air', 'epsilon_medium', 'epsilon_in', 'wavelength_list', 'p', 'FOV_before_windowing');

%% Specify the incident angle or transverse focal location
if strcmpi(incident_profile_type, 'angular')
    incident_angle = -15; % the incident angle. unit: degree
else
    y_f = 0; % transverse focal location. unit: micron
end
%% Construct the permittivity profile
L_air = 50; % thickness of air on the left. unit: micron
W_air = 20; % width of air on the top and bottom. unit: micron
n_left_extra = round(L_air/dx);
n_top_extra = round(W_air/dx);

ny_sca = round(W/dx);
nz_sca = round(L/dx);

nPML = PML.npixels;
n_source = nPML+1; % index of the source plane. The source is placed just in front of the PML.

epsilon = [epsilon_in*ones(nPML*2+ny_sca+2*n_top_extra, n_source+n_left_extra),[epsilon_medium*ones(nPML+n_top_extra, nz_sca);epsilon;epsilon_medium*ones(nPML+n_top_extra, nz_sca)], epsilon_medium*ones(nPML*2+ny_sca+2*n_top_extra, nPML)];
ny_tot = size(epsilon, 1);

%% Specify other data to the syst structure
syst.PML = PML;
syst.dx = dx;
syst.wavelength = 2/(1/wavelength_list(1)+1/wavelength_list(end));
syst.epsilon = epsilon;

%% Construct the input profile
kdx = sqrt(epsilon_in)*2*pi/syst.wavelength*dx;
if strcmpi(incident_profile_type, 'angular')
    % Generate the Planck-taper window that will be applied to the profile.
    % This is used to suppress the source amplitude inside the PML.
    ny_fov = round(FOV_before_windowing/dx);
    win_profile = planckwin(ny_fov, p);

    channels_fov = mesti_build_channels(ny_fov, 'TM', 'periodic', kdx, epsilon_in);
    kydx = channels_fov.kydx_prop;
    % Select the incident angle
    [~, idx] = min(abs(kydx/kdx-sind(incident_angle)));
    E_i_prof = channels_fov.fun_u(kydx(idx)); % the input profile on the focal plane.

    % Generate the list of E^in(z=z_f_air, y).
    E_i = zeros(ny_tot, 1);
    y_i = (1:ny_fov) + round((ny_tot-ny_fov)/2);
    E_i(y_i) = reshape(win_profile, [], 1).*E_i_prof;
else
    % Calculate Gaussian beam wasit for generating profile on the focal plane.
    w0 = syst.wavelength/(pi*NA); % Here, refractive index is already included in NA.

    % Generate the list of E^in(z=z_f_air, y).
    % Here, y.' is an ny-by-1 column vector, and y_f is a 1-by-M_in row vector.
    % So, y.' - y_f is an ny-by-M_in matrix by implicit expansion.
    % Then, E_f is an ny-by-M_in matrix whose m-th column is the cross section
    % of the m-th Gaussian beam at z = z_f_air.
    y = (0.5:ny_tot)*dx-(ny_tot*dx)/2;
    E_i = exp(-(y.' - y_f).^2/(w0^2));
end

% Get properties of propagating channels in the free space for
% projection and backpropagation.
channels = mesti_build_channels(ny_tot, 'TM', 'periodic', kdx, epsilon_in);
u = channels.fun_u(channels.kydx_prop);

% Project E^in(z_f_air, y) onto the propagating channels.
E_i_prop = u'*E_i;

% Backpropagate the profile from z = z_f_air to z = 0 to obtain the
% surface source profile.
z_source = -dx*(n_left_extra+1/2); % the source plane is half a pixel before the interface.
kz = channels.kxdx_prop/dx; % list of wave numbers
E_s_prop = reshape(exp(1i*kz*(z_source-z_f_air)), [], 1).*E_i_prop;

% Determine the line sources.
% In a closed geometry with no PML in y, a line source of
% -2i*nu(a)*u(:,a) generates outgoing waves with transverse profile
% u(:,a). With PML in y, this is not strictly true but is sufficiently
% accurate since E^in(z=0,y) decays exponentially in y.
nu = reshape(channels.sqrt_nu_prop, [], 1).^2;
B_L = u*(nu.*E_s_prop);

%% Specify the input profile
B_struct.data = B_L;
B_struct.pos = [1, n_source, ny_tot, 1];

opts.exclude_PML_in_field_profiles = true;

%% Compute the field profile
% There is no need to use the METIS ordering. The computation time is
% comparable.
%opts.use_METIS = true; 
psi = mesti(syst, B_struct, [], [], opts);

%% Postprocessing
psi = psi(n_top_extra+(1:ny_sca), :); % remove the extra padding on the top and bottom
I = abs(psi).^2; % only save the intensity to reduce the file size
I = I/max(I, [], 'all'); % normalization
I = single(I); % convert the data to single precision

%% Save the field profile data
field_profile_dir = fullfile(data_dir, 'field_profiles');
if ~exist(field_profile_dir, 'dir')
    mkdir(field_profile_dir);
end
save(fullfile(field_profile_dir, ['field_profile_', char(incident_profile_type), '.mat']), 'I', 'L_air');
end
