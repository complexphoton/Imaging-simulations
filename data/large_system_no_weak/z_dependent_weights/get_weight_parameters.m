function weight_params = get_weight_parameters(imaging_method)
% format: ['front_end', 'front_slope', 'front_b', 'back_start', 'back_slope', 'back_b']
% or ['front_end', 'front_slope', 'front_b', 'back_start', 'back_slope', 'back_b', 'peak_pos', 'peak_val']
switch lower(imaging_method)
    case 'smt'
        weight_params = [140; -0.001; 2; 300; 0; 1.8];
    case 'rcm'
        weight_params = [10; -0.0005; 1.5; 300; -0.0005; 0.9];
    case 'oct'
        weight_params = [120; -0.001; 0.5; 300; 0; 0.5];
    case 'ocm'
        weight_params = [120; 0.006; 1.0; 180; -0.004; 1.0; 150; 3.0];
    case 'isam'
        weight_params = [110; 0.005; 2.0; 180; -0.008; 1.9; 150; 3.1];
    case 'synthetic_ocm'
        weight_params = [120; 0.006; 1.0; 180; -0.004; 1.0; 150; 3.0];
    case 'synthetic_oct'
        weight_params = [120; -0.001; 0.5; 300; 0; 0.5];
    otherwise
        error('Unsupported imaging method.')
end

end