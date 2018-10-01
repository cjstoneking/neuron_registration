%% register_single_point 
% For a given point in a 2D slice, finds the best matching point in a 3D
% volume by optimizing local cross-correlation
% Bayesian optimization is used to keep computation time low
% 
%Author: C.J. Stoneking, 2018. cjstoneking@gmail.com

function [output_coords, output_angle, output_error] = register_single_point(ref_volume, probe_volume, ref_coords, probe_radius, varargin)


    %default parameters are based on performance on a test dataset
    %(high-resolution confocal images of brain slices registered to two-photon volume images)
    p = inputParser;
    addParameter(p, 'downsampling_factor', 1);
    single_angle_range = linspace(-pi/4, pi/4, 10);
    addParameter(p, 'angle_range', {single_angle_range, single_angle_range, single_angle_range});
    addParameter(p, 'fractional_cell_offset', 0.25);
    addParameter(p, 'times_to_sample_from_GP', 5);
    addParameter(p, 'antialias', true);
    addParameter(p, 'GP_sig', 1);
    addParameter(p, 'GP_rho', 10000);
    addParameter(p, 'GP_measurement_noise', 0.001);
    addParameter(p, 'proposal_n_init', 10);
    addParameter(p, 'proposal_n_take', 150);

    parse(p, varargin{:});

    %redefine for clarity
    downsampling_factor = p.Results.downsampling_factor;
    angle_range = p.Results.angle_range;
    fractional_cell_offset = p.Results.fractional_cell_offset;
    times_to_sample_from_GP = p.Results.times_to_sample_from_GP;
    antialias = p.Results.antialias;
    GP_sig = p.Results.GP_sig;
    GP_rho = p.Results.GP_rho;
    GP_measurement_noise = p.Results.GP_measurement_noise;
    proposal_n_init = p.Results.proposal_n_init;
    proposal_n_take = p.Results.proposal_n_take;

    for i = 1:3
        angle_step(i) = mean(angle_range{i}(2:end) - angle_range{i}(1:(end-1)));
    end

    probe_volume(isnan(probe_volume)) = 0;
    ref_volume(isnan(ref_volume)) = 0;


    %%
    if(abs(downsampling_factor - 1) > 10^(-6))

        ref_volume = downsample_stack(ref_volume, downsampling_factor);

        probe_volume  = downsample_stack(probe_volume, downsampling_factor);

        ref_coords = ref_coords/downsampling_factor;

        probe_radius = floor(probe_radius/downsampling_factor);

    end

    probe_radius = ceil(probe_radius);

    inner_side_length = floor(sqrt(2)*probe_radius) - 1;
    %side length of each square slice
    %this gives a square which is inscribed
    %in a circle of radius = probe_radius

    rx1 = max(floor(ref_coords(1) - inner_side_length/2), 1);
    rx2 = rx1 + inner_side_length - 1;
    ry1 = max(floor(ref_coords(2) - inner_side_length/2), 1);
    ry2 = ry1 + inner_side_length - 1;

    ref_fft2 = fft2(ref_volume(rx1:rx2, ry1:ry2));

    cell_offset      = 2*fractional_cell_offset*probe_radius;

    grid_coords = get_grid_coords(size(probe_volume), 2*probe_radius, cell_offset, false);
    x0 = grid_coords{1} + probe_radius;
    y0 = grid_coords{2} + probe_radius;
    z0 = grid_coords{3} + probe_radius;

    [x, y, z, a1, a2, a3] = ndgrid(x0, y0, z0, angle_range{1}, angle_range{2}, angle_range{3});

    input_params = [x(:), y(:), z(:), a1(:), a2(:), a3(:)];
    output_params = [];
    error = [];

    %find indices and coordinates of all points in cube that are
    %within circle
    points_in_ball  = zeros(3, (2*probe_radius)^3);
    indices_in_ball = zeros(3, (2*probe_radius)^3);
    pos = 0;

    for i = 1:2*probe_radius
        for j = 1:2*probe_radius
            for k = 1:2*probe_radius
                x = i - probe_radius;
                y = j - probe_radius;
                z = k - probe_radius;
                if((x^2 + y^2 + z^2) <= probe_radius^2)
                    pos = pos + 1;
                    points_in_ball(:, pos)  = [x; y; z];
                    indices_in_ball(:, pos) = [i; j; k];
                end
            end
        end
    end

    points_in_ball = points_in_ball(:,1:pos);
    indices_in_ball = indices_in_ball(:,1:pos);

    %linearize the indices
    indices_in_ball = sub2ind([2*probe_radius, 2*probe_radius, 2*probe_radius],...
        indices_in_ball(1,:), indices_in_ball(2,:), indices_in_ball(3,:));

    if(antialias)
        %compute antialiased rotation by averaging over nearest neighbors
        [off_x, off_y, off_z] = ndgrid(-1:1:1, -1:1:1, -1:1:1);
        offsets = [off_x(:), off_y(:), off_z(:)];
    else
        %no antialias - equivalent to using 1 neighbor = the pixel itself
        offsets = [0, 0, 0];
    end

    indices = {};
    weights = {};
    weights_sum = zeros(1, size(points_in_ball, 2));

    for n = 1:size(offsets, 1)
        indices{n} = zeros(1, size(points_in_ball, 2));
        weights{n} = zeros(1, size(points_in_ball, 2));
    end

    %optimize fft for the current data size
    test_data = randn([inner_side_length, inner_side_length]);
    fftw('dwisdom', []);
    fftw('planner', 'measure');
    fft(test_data);
    fftinfo = fftw('dwisdom');
    fftw('dwisdom', fftinfo);

    i = 0;

    prev_angle = [];

    times_sampled_from_GP = 0;


    while(true)
        i = i + 1;
        if(i > size(input_params, 1))
            if(times_sampled_from_GP >= times_to_sample_from_GP)
                break;
            end
            np = propose_new_params(output_params, error, size(probe_volume), probe_radius, angle_step, times_sampled_from_GP,...
                GP_sig, GP_rho, GP_measurement_noise, proposal_n_init, proposal_n_take);
            if(isempty(np))
                break;
            end
            input_params = [input_params; np];
            times_sampled_from_GP = times_sampled_from_GP + 1;
        end


        angle = input_params(i, 4:6);

        %compute rotation
        %but only if the angle changed from the previous input point
        %for the initial grid, the input points are in blocks of constant angle
        if(isempty(prev_angle) || sum(abs(angle - prev_angle)) > 10^(-6))

            [R1, R2, R3] = rotation_matrices(angle);

            R = R1*R2*R3;
            v = R*points_in_ball + probe_radius;
            %coordinates of the corresponding points
            %in the original image volume

            for n = 1:size(offsets, 1)
                center_pix = floor(v);
                %center pixels in the original image volume
                ind = center_pix + offsets(n,:)';
                ind(ind < 1) = 1;
                ind(ind > 2*probe_radius) = 2*probe_radius;
                lin_ind = sub2ind([2*probe_radius, 2*probe_radius, 2*probe_radius],...
                    ind(1,:), ind(2,:), ind(3,:));
                indices{n} = lin_ind;
                weights{n} = max(sqrt(2*3) - sqrt(sum((v - ind).^2, 1)), 0);
                weights_sum = weights_sum + weights{n};
            end
            prev_angle = angle;
        end

        x1 = floor(input_params(i,1) + 0.5) - probe_radius;
        y1 = floor(input_params(i,2) + 0.5) - probe_radius;
        z1 = floor(input_params(i,3) + 0.5) - probe_radius;
        x2 = x1 + 2*probe_radius - 1;
        y2 = y1 + 2*probe_radius - 1;
        z2 = z1 + 2*probe_radius - 1;

        probe_cube = probe_volume(x1:x2, y1:y2, z1:z2);

        rotated_cube = zeros([2*probe_radius, 2*probe_radius, 2*probe_radius]);
        for n = 1:length(indices)
            rotated_cube(indices_in_ball) = rotated_cube(indices_in_ball)...
                + bsxfun(@times, probe_cube(indices{n}), weights{n}./weights_sum);
        end
        start_ind = ceil((2*probe_radius - inner_side_length)/2);

        rotated_cube_center = rotated_cube(...
            start_ind:(start_ind + inner_side_length - 1),...
            start_ind:(start_ind + inner_side_length - 1),...
            start_ind:(start_ind + inner_side_length - 1));

        lowest_error = inf;
        for z = 1:inner_side_length

            probe_fft2 = fft(fft(rotated_cube_center(:,:,z)).').';

            output = crosscor2D(ref_fft2,probe_fft2);

            if(output(1) <= lowest_error)
                lowest_error = output(1);

                best_local_coords = [inner_side_length/2 - output(2), inner_side_length/2 - output(3)];
                %coords of best matching point within the
                %current rotated square

                current_origin = [x1; y1; z1] + probe_radius;
                %origin of the current cube of the grid

                %r = [-output(3); -output(2); start_ind + z - probe_radius/2];

                %best_probe_coords = (current_origin - R'*r)*downsampling_factor;
                %best_probe_coords = current_origin;
                best_probe_coords_unrotated = [x1, y1, z1] +...
                    [best_local_coords(1), best_local_coords(2), z] + start_ind - [0.5, 0.5, 2];

                r = (best_probe_coords_unrotated' - current_origin);

                best_probe_coords_rotated = current_origin + R*r;
            end
        end

        if(lowest_error < inf)


            output_params(i,:) = [best_probe_coords_rotated(1),...
                best_probe_coords_rotated(2),...
                best_probe_coords_rotated(3),...
                angle(1), angle(2), angle(3),...
                ];
        else
            %never found best coords
            %this happens if ref or probe volume is constant everywhere
            %now the cross-correlation is undefined
            %(because the variance which we normalize by is 0)
            output_params(i,:) = [nan,nan,nan,nan, nan, nan];
        end

        error(i) = lowest_error;
    end

    [best_fit, best_ind] = min(error);

    output_coords = output_params(best_ind, 1:3)*downsampling_factor;
    output_angle  = output_params(best_ind, 4:6);
    output_error  = best_fit;



  end

% register one 2D array to another via crosscorrelation
function [output] = crosscor2D(buf1ft,buf2ft)


    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    error = abs(1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1)));
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,row_shift,col_shift];

end

% get new parameters (ball centers and angles) 
% at which we will evaluate the local cross-correlation during the next step
% of optimization
function [new_params] = propose_new_params(params, error, probe_volume_size, probe_radius, angle_step, search_step_number, GP_sig, GP_rho, GP_measurement_noise, proposal_n_init, proposal_n_take)

    n_initial_points = proposal_n_init;
    %select some number of params with lowest error
    %and perturb these
    
    n_points_to_take = proposal_n_take;
    
    [~, error_sort_ind] = sort(error, 'ascend');
    initial_points = params(error_sort_ind(1:n_initial_points), :);
    %params to initialize search from
    %we perturb these
    
    space_factor = probe_radius/2;
    for i = 1:3
        angle_factor(i) = angle_step(i)/2;
    end
    
    [px, py, pz, pa1, pa2, pa3] = ndgrid((-1:1)*space_factor, (-1:1)*space_factor, (-1:1)*space_factor,...
                           (-1:1)*angle_factor(1), (-1:1)*angle_factor(2), (-1:1)*angle_factor(3));
                       
    perturbations = [px(:), py(:), pz(:), pa1(:), pa2(:), pa3(:)]/(2^search_step_number);

    n_perturbations = size(perturbations, 1);
    
    points_to_predict = zeros(n_initial_points*n_perturbations, 6);
    for p = 1:n_initial_points
        points_to_predict(((p-1)*n_perturbations + 1):(p*n_perturbations),:) = bsxfun(@plus, perturbations, initial_points(p,:));
    end
    within_bounds =  (probe_radius + 1) < points_to_predict(:,1) & points_to_predict(:,1) < (probe_volume_size(1) - probe_radius - 1) &...
                     (probe_radius + 1) < points_to_predict(:,2) & points_to_predict(:,2) < (probe_volume_size(2) - probe_radius - 1) &...
                     (probe_radius + 1) < points_to_predict(:,3) & points_to_predict(:,3) < (probe_volume_size(3) - probe_radius - 1) ;
    
    points_to_predict = points_to_predict(within_bounds, :);         

    predicted_error = GP_predict(params, error, points_to_predict, GP_sig, GP_rho, GP_measurement_noise);
    [~, predicted_error_sort_ind] = sort(predicted_error, 'ascend');

    new_params = points_to_predict(predicted_error_sort_ind(1:min(n_points_to_take, length(predicted_error_sort_ind))),:);
    
end

% Gaussian process (multivariate normal) - based prediction
% as in Rasmussen & Williams 
function [m2] = GP_predict(X1, y1, X2, sig, rho, stability_term)

    if(~iscolumn(y1)); y1=y1';end

    if((size(X2, 1)==1 || size(X2, 2)==1) && ~isrow(X2)); X2 = X2'; end

    added_noise_sd = 0;
    %stability_term = 0.01;
    %sig = 10;
    %rho = 10;

    n1 = size(X1, 1);
    n2 = size(X2, 1);

    %compute distances
    D21 = zeros(n2, n1);
    D11 = zeros(n1, n1);

    for i = 1:n1
       D21(:,i) = sqrt(sum(bsxfun(@minus, X2, X1(i,:)).^2, 2));
       D11(:,i) = sqrt(sum(bsxfun(@minus, X1, X1(i,:)).^2, 2));
    end

    %compute covariances from distances
    C21 = sig^2 * exp(-D21.^2/rho);
    C11 = sig^2 * exp(-D11.^2/rho);

    %compute mean of Gaussian process conditioned on observed values
    m2 = C21*((C11 + (added_noise_sd + stability_term)*eye(n1))\(y1 - mean(y1))) + mean(y1);

end