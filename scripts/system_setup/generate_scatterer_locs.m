function scatterer_loc_list = generate_scatterer_locs(W, L, dx, scatterer_density_list, scatterer_diameter_list, scatterer_min_separation, random_seed_list)
%generate_scatterer_locs Generate the scatterer locations.
%    scatterer_loc_list = generate_scatterer_locs(W, L, dx, scatterer_density_list, scatterer_diameter_list, scatterer_min_separation, random_seed_list)
%
%  === Input Arguments ===
%  W, L (scalar, real): 
%    The width (W) and length (L) of the system in micron.
%  dx (scalar, real):
%    The discretization mesh size in micron.
%  scatterer_density_list (vector, real):
%    The scatterer density in 1/micron^2. Each element corresponds to
%    the density for one type of scatterers.
%  scatterer_diameter_list (vector, real):
%    The scatterer diameter in micron. Each element corresponds to
%    the diameter for one type of scatterers.
%  scatterer_min_separation (symmetric matrix, real):
%    A symmetric matrix M to characterize the scatterer minimum separation.
%    M_{ij} corresponds to the minimum separation between scatterers of
%    type i and j in micron.
%  random_seed_list (vector, integer):
%    The random seed for generating scatterer locations. Each element 
%    corresponds to the random seed for one type of scatterers.
%    
%  === Output Arguments ===
%  scatterer_loc_list (cell array):
%    The scatterer locations. The element scatterer_loc_list{i} is an
%    N-by-2 matrix, where N is the total number of scatterers of type i.
%    Each row of the matrix corresponds to one scatterer location (z, y) in
%    micron.
%
% === Notes ===

if ~issymmetric(scatterer_min_separation)
    error('The scatterer_min_separation must be a symmetric matrix.')
end

n_scatterer_types = length(scatterer_density_list);
scatterer_loc_list = cell(n_scatterer_types, 1);

for i = 1:n_scatterer_types
    rng(random_seed_list(i))
    n_scatterers = round(W*L*scatterer_density_list(i));
    scatterer_locs = zeros(n_scatterers, 2);
    dist_from_boundary = scatterer_diameter_list(i)/2;
    fprintf(['generating type-', num2str(i), ' scatterer locations: ']);
    for j = 1:n_scatterers
        min_distance = zeros(n_scatterer_types, 1); % minimum distance between the generated scatterer and existing scatterers
        touch_boundary = false;
        while any(min_distance < scatterer_min_separation(:, i)) || touch_boundary
            new_scatterer_loc = rand(1, 2).*[L, W];
            % We round the scatterer location here to make the scatterer shape symmetric;
            % The scatterer image could split if its epsilon profile is distorded. 
            new_scatterer_loc = round(new_scatterer_loc/dx)*dx;
            
            % Calculate the minimum distance between the new scatterer and existing scatterers
            for ii = 1:n_scatterer_types
                if ii ~= i
                    scatterer_locs_ii = scatterer_loc_list{ii};
                else
                    scatterer_locs_ii = scatterer_locs(1:j-1, :);
                end

                if ~isempty(scatterer_locs_ii)
                    scatterer_distance = vecnorm(new_scatterer_loc-scatterer_locs_ii, 2, 2);
                    min_distance(ii) = min(scatterer_distance);
                else
                    % Set the distance to inf if scatterers have not
                    % been generated
                    min_distance(ii) = inf;
                end
            end

            if new_scatterer_loc(1) < dist_from_boundary || ...
                    new_scatterer_loc(1) > L - dist_from_boundary || ...
                    new_scatterer_loc(2) < dist_from_boundary || ...
                    new_scatterer_loc(2) > W - dist_from_boundary
                touch_boundary = true;
            else
                touch_boundary = false;
            end
        end
        scatterer_locs(j, :) = new_scatterer_loc;
        textprogressbar(j, j/n_scatterers*100);
    end
    scatterer_loc_list{i} = scatterer_locs;
    fprintf('\n');
end

end