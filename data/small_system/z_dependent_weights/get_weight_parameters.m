function weight_params = get_weight_parameters(imaging_method)
% format: ['front_end', 'front_slope', 'front_b', 'back_start', 'back_slope', 'back_b']
% or ['front_end', 'front_slope', 'front_b', 'back_start', 'back_slope', 'back_b', 'peak_pos', 'peak_val']
switch lower(imaging_method)
    case 'smt'
        weight_params = [60; -0.03; 0.8; 60; 0; 0];
    case 'rcm'
        weight_params = [40; -0.02; 0.6; 60; 0; 0.4];
    case 'oct'
        weight_params = [60; -0.005; 0.5; 60; 0; 0.5];
    case 'ocm'
        weight_params = [20; 0.06; 2; 40; -0.065; 2];
    case 'isam'
        weight_params = [20; 0.01; 2.2; 40; -0.07; 2];
    otherwise
        error('Unsupported imaging method.')
end

end