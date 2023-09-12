clear
data_dir = fullfile('.', 'data', 'large_system');
syst_data_path = fullfile(data_dir, 'system_data.mat');
load(syst_data_path, 'W', 'L');
plot_angular = false;
if plot_angular
    load(fullfile(data_dir, 'field_profiles', 'field_profile_angular.mat'), 'I', 'L_air')
else
    load(fullfile(data_dir, 'field_profiles', 'field_profile_spatial.mat'), 'I', 'L_air')
end

imagesc([-L_air, L], [-W/2, W/2], I)
axis image
colormap(colorcet('L20'));
colorbar
clim([0, 0.4])
xlabel('z (µm)')
ylabel('y (µm)')
set(gca,'TickDir','out');

