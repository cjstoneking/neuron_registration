
function [] = demo_registration()

    % fix seed if needed
    %rng(2020203);

    volume_side_length = 75;
    %side length of the 3D volume

    slice_side_length = 60;
    %side length of the 2D slice

    local_radius = 25;
    %probe radius for registration

    times_sample_from_GP = 5;
    %number of times we sample from Gaussian process

    angle_range_bounds = pi/8;
    angle_discretization_steps = 5;

    fractional_cell_offset = 0.5;

    downsampling_factor = 1;

    GP_sig = 1;
    GP_rho = 10000;
    GP_measurement_noise =  0.001;

    proposal_n_init = 10;
    proposal_n_take = 100;

    VX = generate_test_image(volume_side_length);

    %generate VY via translation + rotation
    VY = nan(volume_side_length)*min(VX(:));
    
    angle = rand(1, 3)*2*angle_range_bounds - angle_range_bounds;
    angle(1) = 0;
    angle(2) = 0;

    [R1, R2, R3] = rotation_matrices(angle);
    R = R1*R2*R3;

    shift = [0, 0, 0];

    [off_x, off_y, off_z] = ndgrid(-1:1:1, -1:1:1, -1:1:1);
    offsets = [off_x(:), off_y(:), off_z(:)];

    for i = 1:size(VY, 1)
        for j = 1:size(VY, 2)
            for k = 1:size(VY, 3)
                x = i - size(VY, 1)/2;
                y = j - size(VY, 2)/2;
                z = k - size(VY, 3)/2;

                v = R'*[x, y, z]' - shift' + [size(VX, 1), size(VX, 2), size(VX, 3)]'/2;
                ind_center = floor(v');

                weights_sum = 0;

                for n = 1:size(offsets, 1)
                    ind = ind_center + offsets(n, :);
                    w   = max(sqrt(3) - sqrt(sum((v' - ind).^2, 2)), 0);
                    weights_sum = weights_sum + w;
                end

                for n = 1:size(offsets, 1)
                    ind = ind_center + offsets(n, :);
                    w   = max(sqrt(3) - sqrt(sum((v' - ind).^2, 2)), 0);

                    if(1 <= ind(1) && ind(1) <= size(VX, 1) &&...
                       1 <= ind(2) && ind(2) <= size(VX, 2) &&...
                       1 <= ind(3) && ind(3) <= size(VX, 3) )
                        if(isnan(VY(i,j,k)))
                            VY(i,j,k) = 0;
                        end
                        VY(i,j,k) = VY(i, j,k) + (w/weights_sum)*VX(ind(1), ind(2), ind(3));
                    end
                end
            end
        end
    end

    VY(isnan(VY)) = 0;

    margin = fix((volume_side_length - slice_side_length)/2);
    VY = VY(margin:( margin + slice_side_length - 1), margin:( margin + slice_side_length - 1));

    y0 = randi(slice_side_length - 2*local_radius, [1 2]) + local_radius;
    a_range = linspace(-angle_range_bounds,...
                        angle_range_bounds,...
                        angle_discretization_steps);


    t0 = tic();

    [x0, ~, ~] = register_single_point(VY, VX,  y0,...
                local_radius,...
                'downsampling_factor', downsampling_factor,...
                'angle_range', {a_range, a_range, a_range},...
                'fractional_cell_offset', fractional_cell_offset,...
                'times_to_sample_from_GP', times_sample_from_GP,...
                'GP_sig', GP_sig,...
                'GP_rho', GP_rho,...
                'GP_measurement_noise', GP_measurement_noise,...
                'proposal_n_init', proposal_n_init,...
                'proposal_n_take', proposal_n_take,...
                'antialias', true);

    computation_time = toc(t0);

    %get the ground-truth value for x0
    ry = [y0(1) - slice_side_length/2; y0(2) - slice_side_length/2; 0];
    x_true = (volume_side_length/2 + R'*ry)';

    euclidean_distance = sqrt(sum((x0 - x_true).^2));
    
    disp(sprintf('computation time:  %f s', computation_time));
    disp(sprintf('distance of estimated point to ground truth: %f pixels', euclidean_distance));
    
    VX_to_plot = VX(margin:( margin + slice_side_length - 1), margin:( margin + slice_side_length - 1), fix(size(VX, 3)/2) ); 
    figure;
    subplot(1,2,1);
    imagesc(VX_to_plot);
    axis equal;
    xlim([0 size(VX_to_plot, 1)]);
    ylim([0 size(VX_to_plot, 2)]);
    title('original image');
    hold on;
    rectangle('Position', [x0(2) - margin - 2, x0(1) - margin - 2, 4, 4],...
                        'Curvature', [1 1], 'EdgeColor', [1 0 0], 'LineWidth', 2);
    subplot(1,2,2);
    imagesc(VY);
    axis equal;
    xlim([0 size(VY, 1)]);
    ylim([0 size(VY, 2)]);
    title('transformed image');
    hold on;
    rectangle('Position', [y0(2)  - 2, y0(1)  - 2, 4, 4],...
                        'Curvature', [1 1], 'EdgeColor', [1 0 0], 'LineWidth', 2);


end



