function weight_params = get_weight_parameters(imaging_method)
% format: ['front_end', 'front_slope', 'front_b', 'back_start', 'back_slope', 'back_b']
% or ['front_end', 'front_slope', 'front_b', 'back_start', 'back_slope', 'back_b', 'peak_pos', 'peak_val']
switch lower(imaging_method)
    case 'smt'
        weight_params = [140; -0.016; 0.2; 300; 0; -1.1];
    case 'rcm'
        weight_params = [50; -0.013; 1.2; 180; -0.0014; 0.55];
    case 'oct'
        weight_params = [120; -0.0072; 0; 240; 0; -0.5];
    case 'ocm'
        weight_params = [100; 0.002; 1.3; 220; -0.002; 0.55];
    case 'isam'
        weight_params = [120; -0.01; 1.1; 180; -0.002; 0.85; 150; 1.50];
    otherwise
        error('Unsupported imaging method.')
end

end