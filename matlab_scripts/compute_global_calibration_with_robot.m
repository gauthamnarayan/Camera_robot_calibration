
function [global_matrix, intrinsics, average_error] = compute_global_calibration_with_robot(robot_poses, depth_image_folder, plane_guess, intrinsics_guess, depth_error_fn)

if (nargin < 4)
    error('usage: [global_matrix, intrinsics, average_error] = compute_global_calibration_with_robot(robot_poses, depth_image_folder, plane_guess, intrinsics_guess, [depth_error_fn])');
elseif (nargin < 5)
    % Found experimentally for one SR300 camera
    depth_error_fn = [0.00034; -0.0002258; 0.00506];
end
 
% rng('default');
% rng(1);

% for i=1:1000
%     w = rand(3,1)*2-1;
%     w_new = rot2skew(skew2rot(w));
%     
%     if (norm(w-w_new)> 0.001)
%         fprintf('[%1.4f %1.4f %1.4f], [%1.4f %1.4f %1.4f]\n', w, w_new);
%     end
% end

global_matrix = zeros(4,3);
intrinsics = zeros(4,1);
average_error = 0;

% We will assume the plane is not perpendicular to the x-y plane, and thus
% the z component is non-zero. In this case we'll assume it's 1.
if (plane_guess(3) == 0)
    % If the user's guess for the plane has it perpendicular to the x-y
    % plane, then just add a bit of a z component
    plane_guess(3) = 0.001;
end

% Make it so the equation of the plane is: a*x + b*y + z + d = 0
plane_guess = plane_guess / plane_guess(3);

% Save robot_poses in case we want to run this code later
save(fullfile(depth_image_folder, 'robot_poses.mat'), 'robot_poses');

% Convert from mm to m
robot_poses(:,1:3) = robot_poses(:,1:3) / 1000.0;
% Make sure quaternions are normalized
robot_poses(:,4:7) = bsxfun(@rdivide, robot_poses(:,4:7), sqrt(sum(robot_poses(:,4:7).^2, 2)));

N = size(robot_poses, 1);

plane_pixels_and_depths = cell(N,1);
plane_equations = cell(N,1);

% Read in depth images, and get a list of pixels that belong to the plane
for i=1:N
    data = load(fullfile(depth_image_folder, sprintf('%d.mat', i-1)));
    depth_image = data.depth;
    [point_cloud, pixels] = deproject_pixel_to_point(depth_image, intrinsics_guess);
    
    ptCloud = pointCloud(point_cloud);
    
    % Now use RANSAC to find a plane relatively close to parallel to the x-y
    % plane (the main wall we're looking at)
    max_distance = 0.005;
    [model1,inlierIndices,outlierIndices] = pcfitplane(ptCloud,...
                max_distance, 'Confidence', 99.9999);
%    length(inlierIndices)
%     figure(1);clf;
%     pcshow(ptCloud.Location(inlierIndices,:), 'b');hold on;
%     pcshow(ptCloud.Location(outlierIndices,:), 'r');drawnow;pause
    
    
    % Take a subset of all the inliers
    idxs = randperm(length(inlierIndices), ceil(length(inlierIndices)/100.0));

    plane = model1.Parameters;
    plane = plane / norm(plane(1:3));
    if (plane(4) < 0)
        plane = -1 * plane;
    end
    
    plane_equations{i} = plane;
    
    % Remember these pixels that belong to the plane
    plane_pixels_and_depths{i} = [pixels(inlierIndices(idxs),:) double(point_cloud(inlierIndices(idxs),3))];
end

% Get a rough guess of the tool frame using plane to plane calibration

P = zeros(N,4);
Q = zeros(N,4);

for i = 1:N
    plane = plane_guess * [quat2rotm(robot_poses(i,4:7)) robot_poses(i,1:3)'; 0 0 0 1];
    plane = plane / norm(plane(1:3));
    if (plane(4) < 0)
        plane = -1 * plane;
    end
    P(i,:) = plane;
    Q(i,:) = plane_equations{i};
end

t = pinv(P(:,1:3)) * (Q(:,4) - P(:,4));
[u,s,v] = svd(P(:,1:3)' * Q(:,1:3));
R = u*v';

transform_guess = [t' rotm2quat(R)];

fprintf('tool frame guess: %1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f\n', transform_guess);



% 
% for i = 1:N
%     plane = plane_guess * [quat2rotm(robot_poses(i,4:7)) robot_poses(i,1:3)'; 0 0 0 1] * [R t;0 0 0 1];
%     plane = plane / norm(plane(1:3));
%     if (plane(4) < 0)
%         plane = -1 * plane;
%     end
%     %fprintf('%1.4f %1.4f %1.4f %1.4f, %1.4f %1.4f %1.4f %1.4f\n', plane, plane_equations{i});
%     fprintf('%1.4f %1.4f %1.4f %1.4f\n', plane_equations{i} * inv([quat2rotm(robot_poses(i,4:7)) robot_poses(i,1:3)'; 0 0 0 1] * [R t;0 0 0 1]));
% end


% Setup optimization

x0 = [0.00187659, 1.02581, -0.107771, 0.00463381, 0.9949, -0.127648, 0.00426012, 1.00811, -0.113047, 1,1,0,0, plane_guess([1 2 4]), t', rot2skew(R)']';
%x0 = [0, 1, 0, 0, 1, 0, 0, 1, 0, 1,1,0,0, plane_guess([1 2 4]), t', rot2skew(R)']';

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 100000, 'MaxIterations', 4000, 'Display', 'iter');
[params, resnorm] = lsqnonlin(@(x)global_match_residuals(x(1:9), x(10:13), x(14:16), x(17:22), robot_poses, plane_pixels_and_depths, intrinsics_guess', depth_error_fn), x0, [], [], options);

average_error = resnorm / N;

% Return results
global_matrix_partial = params(1:9);
global_matrix = [global_matrix_partial(1) global_matrix_partial(2) global_matrix_partial(3);
    global_matrix_partial(4) global_matrix_partial(5) global_matrix_partial(6);
    global_matrix_partial(7) global_matrix_partial(8) global_matrix_partial(9);
    0 0 0];
global_matrix(4,:) = global_matrix(2,:) + global_matrix(3,:) - global_matrix(1,:);

intrinsics = zeros(4,1);
intrinsics(1) = intrinsics_guess(1) * params(10);
intrinsics(2) = intrinsics_guess(2) * params(11);
intrinsics(3) = intrinsics_guess(3) + params(12);
intrinsics(4) = intrinsics_guess(4) + params(13);

base_plane_params = params(14:16);
base_plane = [base_plane_params(1) base_plane_params(2) 1.0 base_plane_params(3)];

tool_transform = params(17:22);
tool_matrix = [skew2rot(tool_transform(4:6)) [tool_transform(1); tool_transform(2); tool_transform(3)]; 0 0 0 1];

fprintf('fx: %1.4f, fy: %1.4f, cx: %1.4f, cy: %1.4f\n', intrinsics);
fprintf('base plane: [%1.4f, %1.4f, %1.4f, %1.4f]\n', base_plane);
fprintf('tool transform: [%1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f]\n', tool_transform(1:3), rotm2quat(tool_matrix(1:3, 1:3)));
fprintf('global matrix:\n');
for i=1:4
    fprintf('%1.6f %1.6f %1.6f\n', global_matrix(i,:));
end


% x0 = [1,1,0,0, plane_guess([1 2 4]), t', rot2skew(R)']';
% options = optimset('MaxFunEvals', 2000*13, 'MaxIter', 2000*13, 'Display', 'iter');
% params = fminunc(@(x)global_match_cost(x(1:4), x(5:7), x(8:13), robot_poses, plane_pixels_and_depths, intrinsics_guess', depth_error_fn), x0, options);

% x0 = params(10:22);
% options = optimset('MaxFunEvals', 2000*13, 'MaxIter', 2000*13, 'Display', 'iter');
% params = fminsearch(@(x)global_match_cost(x(1:4), x(5:7), x(8:13), robot_poses, plane_pixels_and_depths, intrinsics_guess', depth_error_fn), x0, options);
% 
% % Return results
% [cost, global_matrix] = global_match_cost_with_global_matrix(params(1:4), params(5:7), params(8:13), robot_poses, plane_pixels_and_depths, intrinsics_guess', depth_error_fn);
% 
% intrinsics(1) = intrinsics_guess(1) * params(1);
% intrinsics(2) = intrinsics_guess(2) * params(2);
% intrinsics(3) = intrinsics_guess(3) + params(3);
% intrinsics(4) = intrinsics_guess(4) + params(4);
% 
% base_plane_params = params(5:7);
% base_plane = [base_plane_params(1) base_plane_params(2) 1.0 base_plane_params(3)];
% 
% tool_transform = params(8:13);
% tool_matrix = [skew2rot(tool_transform(4:6)) [tool_transform(1); tool_transform(2); tool_transform(3)]; 0 0 0 1];
% 
% fprintf('fx: %1.4f, fy: %1.4f, cx: %1.4f, cy: %1.4f\n', intrinsics);
% fprintf('base plane: [%1.4f, %1.4f, %1.4f, %1.4f]\n', base_plane);
% fprintf('tool transform: [%1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f]\n', tool_transform(1:3), rotm2quat(tool_matrix(1:3, 1:3)));
% fprintf('global matrix:\n');
% for i=1:4
%     fprintf('%1.6f %1.6f %1.6f\n', global_matrix(i,:));
% end

end

function cost = global_match_cost(intrinsics_delta, base_plane_params, tool_transform, robot_poses, plane_pixels_and_depths, intrinsics, depth_error_fn)

cost = global_match_cost_with_global_matrix(intrinsics_delta, base_plane_params, tool_transform, robot_poses, plane_pixels_and_depths, intrinsics, depth_error_fn);

end

function [cost, global_matrix] = global_match_cost_with_global_matrix(intrinsics_delta, base_plane_params, tool_transform, robot_poses, plane_pixels_and_depths, intrinsics, depth_error_fn)

% Step 1: Extract things from arguments

% Number of camera images
N = size(robot_poses,1);

% Tool transform from [x y z wx wy wz]
tool_matrix = [skew2rot(tool_transform(4:6)) [tool_transform(1); tool_transform(2); tool_transform(3)]; 0 0 0 1];

% Camera intrinsics
fx = intrinsics(1) * intrinsics_delta(1);
fy = intrinsics(2) * intrinsics_delta(2);
cx = intrinsics(3) + intrinsics_delta(3);
cy = intrinsics(4) + intrinsics_delta(4);

% Base plane
base_plane = [base_plane_params(1) base_plane_params(2) 1.0 base_plane_params(3)];

% Each column of this matrix is the (u,v) coordinate corresponding to each
% row of the global matrix
st = [1 641 1 641;
      1 1 481 481];

A = [];
b = [];
w = [];
  
for i=1:N
    % Find equation of plane in camera frame
    plane_in_camera = base_plane * [quat2rotsingle(robot_poses(i,4:7)) robot_poses(i,1:3)'; 0 0 0 1] * tool_matrix;
    
    % Correct depths according to our global matrix for each of our pixels
    % of interest
    u = plane_pixels_and_depths{i}(:,1);
    v = plane_pixels_and_depths{i}(:,2);
    z = plane_pixels_and_depths{i}(:,3);
    
    % Get weights and function values at each of the corners
    w_u = 1 - abs(bsxfun(@minus, u, st(1,:))) / 640;
    w_v = 1 - abs(bsxfun(@minus, v, st(2,:))) / 480;
    
    A = [A; bsxfun(@times, (w_u(:,1) .* w_v(:,1) - w_u(:,4) .* w_v(:,4)), [ones(length(z),1), z, z.^2]), bsxfun(@times, (w_u(:,2) .* w_v(:,2) + w_u(:,4) .* w_v(:,4)), [ones(length(z),1), z, z.^2]), bsxfun(@times, (w_u(:,3) .* w_v(:,3) + w_u(:,4) .* w_v(:,4)), [ones(length(z),1), z, z.^2])];

    % Now, using our intrinsics, find the distance for each of our pixels
    % of interest to the plane
    depths = (-plane_in_camera(4) ./ (plane_in_camera(1:3) * [(u'-cx)/fx; (v'-cy)/fy; ones(1,length(u))]))';
    
    b = [b; depths];

    % The difference between the distance to the plane and our corrected
    % depth will be our cost function
    depth_error = [ones(length(z),1), z, z.^2] * depth_error_fn;
    w = [w; 1 ./ (length(z) * depth_error.^2)];
end

x = lscov(A,b,w);

cost = sum(w .* (A*x - b).^2);

if(nargout == 2)
    global_matrix = [x(1) x(2) x(3);x(4) x(5) x(6); x(7) x(8) x(9); 0 0 0];
    global_matrix(4,:) = global_matrix(2,:) + global_matrix(3,:) - global_matrix(1,:);
end
    

end

function residuals = global_match_residuals(global_matrix_partial, intrinsics_delta, base_plane_params, tool_transform, robot_poses, plane_pixels_and_depths, intrinsics, depth_error_fn)

% Step 1: Extract things from arguments

% Number of camera images
N = size(robot_poses,1);

% Tool transform from [x y z wx wy wz]
tool_matrix = [skew2rot(tool_transform(4:6)) [tool_transform(1); tool_transform(2); tool_transform(3)]; 0 0 0 1];

% Camera intrinsics
fx = intrinsics(1) * intrinsics_delta(1);
fy = intrinsics(2) * intrinsics_delta(2);
cx = intrinsics(3) + intrinsics_delta(3);
cy = intrinsics(4) + intrinsics_delta(4);

% Base plane
base_plane = [base_plane_params(1) base_plane_params(2) 1.0 base_plane_params(3)];

% Global polynomial matrix, where each row is [a b c], a + b * z + c * z^2
global_matrix = [global_matrix_partial(1) global_matrix_partial(2) global_matrix_partial(3);
    global_matrix_partial(4) global_matrix_partial(5) global_matrix_partial(6);
    global_matrix_partial(7) global_matrix_partial(8) global_matrix_partial(9);
    0 0 0];
global_matrix(4,:) = global_matrix(2,:) + global_matrix(3,:) - global_matrix(1,:);


% Each column of this matrix is the (u,v) coordinate corresponding to each
% row of the global matrix
st = [1 641 1 641;
      1 1 481 481];

residuals = [];
  
for i=1:N
    % Find equation of plane in camera frame
    plane_in_camera = base_plane * [quat2rotsingle(robot_poses(i,4:7)) robot_poses(i,1:3)'; 0 0 0 1] * tool_matrix;
    
    % Correct depths according to our global matrix for each of our pixels
    % of interest
    u = plane_pixels_and_depths{i}(:,1);
    v = plane_pixels_and_depths{i}(:,2);
    z = plane_pixels_and_depths{i}(:,3);
    
    % Get weights and function values at each of the corners
    w_u = 1 - abs(bsxfun(@minus, u, st(1,:))) / 640;
    w_v = 1 - abs(bsxfun(@minus, v, st(2,:))) / 480;
    
    g = [ones(length(z),1), z, z.^2] * global_matrix';
    
    % Correct our depth via a weighted sum
    new_z = sum(w_u .* w_v .* g,2);
    
    % Now, using our intrinsics, find the distance for each of our pixels
    % of interest to the plane
    depths = (-plane_in_camera(4) ./ (plane_in_camera(1:3) * [(u'-cx)/fx; (v'-cy)/fy; ones(1,length(u))]))';
    
    % The difference between the distance to the plane and our corrected
    % depth will be our cost function
    depth_error = [ones(length(z),1), z, z.^2] * depth_error_fn;
    residuals = [residuals; (new_z - depths) ./ (sqrt(length(z))*depth_error)];
end


end

% Convert a [wx;wy;wz] which are the coordinates of a skew symmetric matrix
% into a rotation matrix
function R = skew2rot(w)
    nw = norm(w);
    if (norm(w) == 0)
        R = eye(3);
    else
        w_hat = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
        R = eye(3) + sin(nw)/nw * w_hat + (1 - cos(nw))/nw^2 * (w_hat * w_hat);
    end
end

% Convert a Rotation matrix into skew symmetric coordinates. Note that this
% actually returns w*theta, where w is a unit vector. 
function w = rot2skew(R)

theta = acos((trace(R) - 1) / 2);
if (theta == 0)
    w = [0;0;0];
else
    w = 1/(2*sin(theta)) * [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)] * theta;
end
end

function R = quat2rotsingle(q)
R = [q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4), 2*(q(1)*q(3) - q(1)*q(4)),  2*(q(2)*q(4) + q(1)*q(3));
    2*(q(1)*q(3) + q(1)*q(4)), q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4), 2*(q(3)*q(4) - q(1)*q(2));
    2*(q(2)*q(4) - q(1)*q(3)), 2*(q(3)*q(4) + q(1)*q(2)), q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)];
end


% Given a depth image and intrinsics, return a point cloud
% depth_image: a h x w matrix(image) of depth values in meters
% intrinsics: intrinsics of the depth camera: [fx;fy;cx;cy]
% point_cloud: a (h x w) x 3 matrix of a point cloud, where the first
% dimension is row-major order, and the second dimension is [x,y,z]
function [point_cloud, pixels] = deproject_pixel_to_point(depth_image, intrinsics)

fx = intrinsics(1);
fy = intrinsics(2);
cx = intrinsics(3);
cy = intrinsics(4);

% Get the pixel values for depth image indices
[pixels_x, pixels_y] = meshgrid(1:size(depth_image,2), 1:size(depth_image,1));

% Deproject pixel (note that we're assuming no distortion)
px = (pixels_x - cx) / fx .* depth_image;
py = (pixels_y - cy) / fy .* depth_image;
pz = depth_image;

% Now form a (h x w) x 3 matrix which represents our point cloud
point_cloud = [px(:) py(:) pz(:)];
pixels = [pixels_x(:), pixels_y(:)];

bad_idx = (pz(:) < 0.001);
point_cloud(bad_idx,:) = [];
pixels(bad_idx, :) = [];

% point_cloud = zeros(size(depth_image,1), size(depth_image,2), 3);
% point_cloud(:,:,1) = px;
% point_cloud(:,:,2) = py;
% point_cloud(:,:,3) = pz;

end