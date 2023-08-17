function scatterer_loc_list = generate_scatterer_locs(W, L, dx, scatterer_density_list, scatterer_diameter_list, scatterer_min_separation, random_seed_list)
% scatterer_loc_list: a length-N cell array.
% scatterer_density_list: a length-N column vector;
% scatterer_diameter_list: a length-N column vector;
% scatterer_min_separation: an N-by-N symmetric matrix.

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
    fprintf(['generating scatterer locations for type ', num2str(i), ' : ']);
    for j = 1:n_scatterers
        min_distance = zeros(n_scatterer_types, 1); % minimum distance between the generated scatterer and existing scatterers
        touch_boundary = 0;
        while any(min_distance < scatterer_min_separation(:, i)) || touch_boundary
            new_scatterer_loc = rand(1, 2).*[L, W];
            % We round the scatterer location here to make the scatterer shape symmetric;
            % The scatterer image could split if its epsilon profile is distorded. 
            new_scatterer_loc = round(new_scatterer_loc/dx)*dx;
            
            %% Calculate the minimum distance between the new scatterer and existing scatterers
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
                touch_boundary = 1;
            else
                touch_boundary = 0;
            end
        end
        scatterer_locs(j, :) = new_scatterer_loc;
        textprogressbar(j, j/n_scatterers*100);
    end
    scatterer_loc_list{i} = scatterer_locs;
    fprintf('\n');
end

end