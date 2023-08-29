data_dir = fullfile('.', 'data', 'large_system');
%data_dir = fullfile('.', 'data', 'large_system_no_weak');
tic
recon_smt(data_dir); % reconstruct the SMT image
recon_rcm(data_dir); % reconstruct the RCM image
recon_ocm_or_oct(data_dir, 'ocm'); % reconstruct the OCM image
recon_ocm_or_oct(data_dir, 'oct'); % reconstruct the OCT image
recon_isam(data_dir); % reconstruct the ISAM image
toc