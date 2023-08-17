clear
%% Set the system for plotting
% 0: the default small system
% 1: the large system from the paper
% 2: the large system without low-index-contrast scatterers
system_idx = 1;

save_figs = false;
combine_figs = false;
show_titles_and_labels = true;
show_colorbars = false;

apply_weight = true;
overlap_target_ground_truth = false;
highlight_ref = false;
draw_boxes = false;

font_size = 8;
font_name = 'Arial';

switch system_idx
    case 0
        data_dir = fullfile('.', 'data', 'small_system');
        figs_dir = fullfile('.', 'figs', 'small_system', 'reconstructed_images');
    case 1
        data_dir = fullfile('.', 'data', 'large_system');
        figs_dir = fullfile('.', 'figs', 'large_system', 'reconstructed_images');
    case 2
        data_dir = fullfile('.', 'data', 'large_system_no_weak');
        figs_dir = fullfile('.', 'figs', 'large_system_no_weak', 'reconstructed_images');
    otherwise
        error("Undefined system index!")
end

syst_params_filename = fullfile(data_dir, 'system_data.mat');
load(syst_params_filename, 'FOV', 'W', 'L', 'W_image', 'L_image', 'epsilon', 'dx', 'z_f_bg', 'NA', ...
    'dy_image', 'dz_image', 'target_locs', 'epsilon_medium', 'PML');

plotting_list = {'Ground_truth', 'SMT', 'RCM', 'OCT', 'OCM', 'ISAM'};
%plotting_list = {'OCT', 'Synthetic_OCT', 'OCM', 'Synthetic_OCM'};
%plotting_list = {'ISAM'};
%plotting_list = {'Ground_truth'};

% colorbar range. xxx_idx: colorbar for the system with index=idx
clim_dict = dictionary('ground_truth_0', {[0, 1]}, 'smt_0', {[0, 3]}, ...
    'rcm_0', {[0, 2]}, 'oct_0', {[0, 4]}, 'ocm_0', {[0, 4]}, 'isam_0', {[0, 5]}, ...
    'ground_truth_1', {[0, 1]}, 'smt_1', {[0, 1.5]}, 'rcm_1', {[0, 1.2]}, ...
    'oct_1', {[0, 2]}, 'ocm_1', {[0, 1.4]}, 'isam_1', {[0, 0.9]}, ...
    'synthetic_oct_1', {[0, 1]}, 'synthetic_ocm_1', {[0, 1.4]}, ...
    'ground_truth_2', {[0, 1]}, 'smt_2', {[0, 1]}, 'rcm_2', {[0, 1.2]}, ...
    'oct_2', {[0, 2]}, 'ocm_2', {[0, 1.0]}, 'isam_2', {[0, 1.2]}, ...
    'synthetic_oct_2', {[0, 2]}, 'synthetic_ocm_2', {[0, 0.8]});

if overlap_target_ground_truth
    filter_scatterer_locs = @(locs, z_range, y_range) z_range(1) < locs(:, 1) & ...
        locs(:, 1) < z_range(end) & y_range(1) < locs(:, 2) & locs(:, 2) < y_range(end);
    target_locs_idx = filter_scatterer_locs(target_locs, [0, L_image], [-W_image/2, W_image/2]+W/2);
    target_locs = target_locs(target_locs_idx, :);
end

if combine_figs
    tiledlayout(3, 2, TileSpacing="compact")
    set(gcf,'Renderer', 'painters', 'color', 'w', 'Position', [10 10 500 750])
end

% Generate the ground truth
if any(strcmpi(plotting_list, 'ground_truth'))
    fprintf('generating the ground truth...\n');
    n_scatterers = length(target_locs);
    epsilon_bg = 0.0; epsilon_scatterer = 1.0^2*ones(n_scatterers, 1); target_diam = 1.0*ones(n_scatterers, 1);
    epsilon = generate_permittivity_profile(W, L, dy_image, dz_image, epsilon_bg, target_locs, target_diam, epsilon_scatterer);
end

for ii = 1:length(plotting_list)

    if combine_figs
        nexttile
    else
        figure(ii)
    end

    plotting_name = plotting_list{ii};
    if ~strcmpi(plotting_name, 'ground_truth')
        if apply_weight
            load(fullfile(data_dir, 'z_dependent_weights', lower(plotting_name)), 'weight');
        else
            weight = 1;
        end

        load(fullfile(data_dir, 'reconstructed_images', lower(plotting_name)), 'y_image', 'z_image', 'I');
        imagesc(z_image, y_image, flipud(I./weight));
    else
        %% plot the permittivity profile of high-index scatterers
        nz = round(L/dz_image); ny = round(W/dy_image);
        z_epsilon = (0.5:nz)*dz_image;
        y_epsilon = (0.5:ny)*dy_image;
        imagesc(z_epsilon, y_epsilon-W/2, flipud(epsilon))

        if draw_boxes
            rectangle('Position', [50, -75, 50, 150], 'LineWidth', 0.5, 'EdgeColor', '#ee7733');
            rectangle('Position', [125, -75, 50, 150], 'LineWidth', 0.5, 'EdgeColor', '#009988');
            rectangle('Position', [200, -75, 50, 150], 'LineWidth', 0.5, 'EdgeColor', '#ee3377');
        end
    end

    axis image
    colormap(colorcet_L3());
    clim(cell2mat(clim_dict([lower(plotting_name), '_', num2str(system_idx)])));

    set(gca,'TickDir','out');
    set(gca, 'YDir','normal');
    xlim([0, L_image]);
    ylim([-W_image/2, W_image/2]);
    xticks(0:L_image/3:L_image);
    yticks(-W_image/2:W_image/3:W_image/2);

    if overlap_target_ground_truth
        hold on
        % flip the y coordinate of target_locs here
        plot(target_locs(:, 1), -(target_locs(:, 2)-W/2), 'bo', 'MarkerSize', 5, 'LineWidth', 2);
    end

    if highlight_ref
        hold on
        plot(z_f_bg*[1, 1] , W_image/2*[-1, 1], 'w--');
    end

    if show_titles_and_labels
        title(strrep(plotting_name, '_', ' '));
        xlabel('{\it z} (µm)');
        ylabel('{\it y} (µm)');
        set(gca, 'fontname', font_name, 'fontsize', font_size);
    else
        set(gca, 'XTickLabel', [], 'YTickLabel', [])
    end

    if show_colorbars
        c = colorbar;
        c.TickLength = 0;
        c.Ticks = [];
    end

    if ~combine_figs && save_figs
        if ~exist(figs_dir, 'dir')
            mkdir(figs_dir)
        end
        exportgraphics(gcf, fullfile(figs_dir, [lower(plotting_name), '.jpg']), 'Resolution', 300)
    end
end

if combine_figs && save_figs
    if ~exist(figs_dir, 'dir')
        mkdir(figs_dir)
    end
    exportgraphics(gcf, fullfile(figs_dir, 'combined_results.jpg'), 'Resolution', 300)
end
