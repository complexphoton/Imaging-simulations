function epsilon = generate_permittivity_profile(W, L, dy, dz, epsilon_medium, scatterer_locs, scatterer_diam, scatterer_epsilon)
ny = W/dy; nz = L/dz;
subpixel_size = 5;

nz_temp = nz*subpixel_size;
ny_temp = ny*subpixel_size;

n_scatterers = length(scatterer_locs);

dy_temp = dy/subpixel_size;
dz_temp = dz/subpixel_size;
epsilon_temp = epsilon_medium*ones(ny_temp,nz_temp);
fprintf('generating epsilon on a finer grid: ')
for ns = 1:n_scatterers
    textprogressbar(ns, ns/n_scatterers*100)
    z0 = scatterer_locs(ns, 1);
    y0 = scatterer_locs(ns, 2);
    rad0 = scatterer_diam(ns)/2;

    % location of (x0,y0) in terms of index (n from 1 to nz, m from 1 to ny)
    n0 = 0.5 + z0/dz_temp;
    m0 = 0.5 + y0/dy_temp;
    % create a box surrounding (x0,y0) with edge length = 2*rad0
    dn = rad0/dz_temp; dm = rad0/dy_temp;
    n1 = max([1,  round(n0 - dn)]);
    n2 = min([nz_temp, round(n0 + dn)]);
    m1 = max([1,  round(m0 - dm)]);
    m2 = min([ny_temp, round(m0 + dm)]);
    z_local = ((n1:n2)-0.5)*dz_temp;
    y_local = ((m1:m2)-0.5)*dy_temp;
    [Z_local,Y_local] = meshgrid(z_local,y_local);
    % fill in the scatterer within this local box
    % this step would be very slow if we use the full [Z,Y] instead of [Z_local,Y_local]
    [m_local, n_local] = find((((Z_local-z0).^2 + (Y_local-y0).^2)) <= (rad0^2));
    % conver to linear index in the global coordinate
    ind_n = sub2ind([ny_temp, nz_temp], m_local+(m1-1), n_local+(n1-1));
    epsilon_temp(ind_n) = scatterer_epsilon(ns);
end
fprintf('\n')

fprintf('subpixel smoothing: ')
%% Average inside each pixel
% average over y
epsilon_avg_over_y = zeros(ny,nz_temp);
for n = 1:nz_temp
    textprogressbar(n, n/(nz_temp+ny)*100)
    epsilon_avg_over_y(:,n) = mean(reshape(epsilon_temp(:,n), subpixel_size, ny), 1).';
end
% average over z
epsilon = zeros(ny,nz);
for m = 1:ny
    textprogressbar(m+nz_temp, (m+nz_temp)/(nz_temp+ny)*100)
    epsilon(m,:) = mean(reshape(epsilon_avg_over_y(m,:).', subpixel_size, nz), 1);
end
fprintf('done\n')
end