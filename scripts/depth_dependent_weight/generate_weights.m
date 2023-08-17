clear

data_dir = fullfile('.', 'data', 'small_system');
data_dir = fullfile('.', 'data', 'large_system');
%data_dir = fullfile('.', 'data', 'large_system_no_weak');

% specify imaging method: 'smt', 'rcm', 'oct', 'ocm', 'isam'
imaging_method = 'smt';

plot_figs = true;
save_weight = false;

load(fullfile(data_dir, 'reconstructed_images', imaging_method), 'y_image', 'z_image', 'I');
weight_dir = fullfile(data_dir, 'z_dependent_weights');

if ~isfile(fullfile(weight_dir, 'get_weight_parameters.m'))
    error(['Cannot find get_weight_parameters.m here: ', weight_dir,'! Please define the weight parameters first.'])
end

curr_dir = cd(weight_dir);
weight_params = get_weight_parameters(imaging_method);
cd(curr_dir);

% Obtain max image intensity at each depth
max_intensity_at_each_depth = max(I);

if length(weight_params) == 6
        % middle part: 3-rd order polynomial.
    M = [weight_params(1).^(3:-1:0);...
        [3*weight_params(1)^2, 2*weight_params(1), 1, 0];...
        weight_params(4).^(3:-1:0);
        [3*weight_params(4)^2, 2*weight_params(4), 1, 0]];
    RHS = [weight_params(3); weight_params(2); weight_params(6); weight_params(5)];
elseif length(weight_params) == 8
        % Fit the peak with a 5-th order polynomial.
    M = [weight_params(1).^(5:-1:0);...
        [weight_params(1).^(4:-1:0), 0].*(5:-1:0);...
        weight_params(4).^(5:-1:0);...
        [weight_params(4).^(4:-1:0), 0].*(5:-1:0);...
        weight_params(7).^(5:-1:0);...
        [weight_params(7).^(4:-1:0), 0].*(5:-1:0)];
    RHS = [weight_params(3);weight_params(2);weight_params(6);weight_params(5);weight_params(8);0];
else
    error('Invalid weight parameters!')
end
coefs = pinv(M)*RHS;

%% Generate the weight
weight = zeros(1, length(z_image));

front_range = z_image(z_image<=weight_params(1));
middle_range = z_image(z_image>weight_params(1)&z_image<=weight_params(4));
back_range = z_image(z_image>weight_params(4));

weight(z_image<=weight_params(1)) = 10.^(weight_params(2)*(front_range-weight_params(1))+weight_params(3));
weight(z_image>weight_params(4)) = 10.^(weight_params(5)*(back_range-weight_params(4))+weight_params(6));
weight(z_image>weight_params(1)&z_image<=weight_params(4)) = 10.^polyval(coefs, middle_range);


if plot_figs
    semilogy(z_image, max_intensity_at_each_depth, 'k');
    hold on
    semilogy(z_image, weight, 'b', LineWidth=2);

    legend('I_{\rm max} at each depth', 'weight')
    xlabel('z (\mum)')
    set(gca, 'fontsize', 16);
    ylim([1e-3, 1e4])
end

if save_weight
    save(fullfile(weight_dir, imaging_method), 'weight', 'max_intensity_at_each_depth', 'z_image');
end