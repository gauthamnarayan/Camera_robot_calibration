%
% rs_extrinsic_cal.m
%
% Author: Robbie Paolini, CMU Robotics Institute
%
% Last Modified: 11/21/2017
%
% This file finds camera poses of an arbitrary number of depth cameras that 
% are pointed at a workspace that a robot moves in. The robot holds a sphere
% of known radius at a given distance away from its tool frame. The robot
% moves to various poses in the workspace, and each of the cameras save one
% point cloud file at each pose. By finding the center of the sphere in
% these point clouds, we can correlate the pose of the camera with the base
% frame of the robot. 
%
% This file takes as input a folder generated from extrinsic_cal.cpp from
% the rs_extrinsic_cal ROS package. The inside of the passed folder should
% be as follows:
% * a robot_poses.txt file, which is a N x 7 matrix of [x y z q0 qx qy qz]
% robot poses with identity work and tool frames.
% M folders with titles of camera names, and within each folder, N .pcd
% files from 1.pcd to N.pcd, which contain point cloud files taken by each
% camera when the robot was at each pose.
%
% The output of the program is a N length struct with the following fields:
% camera_name: string name of camera
% camera_frame: [x y z q0 qx qy qz] pose of camera with respect to robot
% base
% marker_position: [x y z] center of sphere with respect to robot tool
% frame
% average_error: The average error of the sphere center in each point cloud
% image. This number should be as small as possible.
%
%

function [robot_poses, camera_positions, trimmed_clouds] = rs_extrinsic_cal(folder)
%function [camera_frames, marker_position, average_errors] = rs_extrinsic_cal(folder)

% Assuming current MLab directory structure, add the package
% extrinsic_cal's matlab folder to matlab's path
path = mfilename('fullpath');
idxs = strfind(path, '/');
extrinsic_cal_path = fullfile(path(1:idxs(end-3)), 'extrinsic_cal', 'matlab_scripts');
addpath(extrinsic_cal_path);


% Known radius of sphere in units of point clouds and robot_poses.txt
SPHERE_RAD = 0.02853; % Pool cue. % 0.075/2 wooden ball

% The number of pixels per squared radians. Used for computing how many
% points on a sphere we should see if it is a given distance away. This has
% to do with the field of view and pixel resolution of the depth camera
PIX_PER_SQRAD = 27000 / (91.2/180*pi) / (65.0/180*pi);

% The minimum angle of the sphere we'll allow a fit for
VIEW_ANGLE = pi/6;

% Rough calibration of each camera using voxel background subtraction, 
% ICP, and a forgiving RANSAC
[camera_calibrations, all_clouds, robot_poses] = rough_cal(folder, SPHERE_RAD, VIEW_ANGLE, PIX_PER_SQRAD);


fprintf('after rough calibration:\n');
for i=1:length(camera_calibrations)
    fprintf('%1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f\n', camera_calibrations{i}.camera_frame);
end
fprintf('\n');

% Fine calibration of each camera using the rough estimate to do better
% background subtraction, ICP, and a less forgiving RANSAC
[camera_calibrations, trimmed_clouds] = fine_cal(camera_calibrations, all_clouds, robot_poses, SPHERE_RAD, VIEW_ANGLE, PIX_PER_SQRAD);

fprintf('after fine calibration:\n');
for i=1:length(camera_calibrations)
    fprintf('%1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f\n', camera_calibrations{i}.camera_frame);
    pose = camera_calibrations{i}.camera_frame;
    pose([4 5 6 7]) = pose([5 6 7 4]);
    fprintf('For ROS: %1.4f %1.4f %1.4f %1.5f %1.5f %1.5f %1.5f\n', pose);
end
fprintf('\n');

robot_poses = cell(length(camera_calibrations),1);
camera_positions = cell(length(camera_calibrations),1);
for i=1:length(robot_poses)
    robot_poses{i} = camera_calibrations{i}.robot_poses;
    camera_positions{i} = camera_calibrations{i}.camera_positions;
end

end

function [camera_frames, marker_position, average_errors] = radius_cal(robot_poses, camera_positions)

options = optimoptions('fminunc','MaxIterations', 4000, 'MaxFunctionEvaluations', 3000, 'Display', 'off');
dr = fminunc(@(dr)radius_cost(robot_poses, camera_positions,dr),0,options);

camera_positions = correct_camera_positions(camera_positions, dr);

[camera_frames, marker_position, average_errors] = multi_extrinsic_cal(robot_poses, camera_positions);

end

function cost = radius_cost(robot_poses, camera_positions, dr)

camera_positions = correct_camera_positions(camera_positions, dr);
[camera_frames, marker_position, average_errors] = multi_extrinsic_cal(robot_poses, camera_positions);
cost = mean(average_errors);

end

function new_camera_positions = correct_camera_positions(camera_positions, dr)
N = length(camera_positions);
new_camera_positions = cell(N,1);

for i=1:N
    % compute distance of each point from camera center
    dists = sqrt(sum(camera_positions{i}.^2,2));
    
    % Now simply change the length of that vector by the requested
    % correction amount
    new_camera_positions{i} = bsxfun(@times, bsxfun(@rdivide, camera_positions{i}, dists), (dists + dr)); 
end

end

% function [camera_calibrations] = global_cal(camera_calibrations, trimmed_clouds, robot_poses, SPHERE_RAD)
% 
% NUM_CAMERAS = size(trimmed_clouds,1);
% x0 = zeros(NUM_CAMERAS*7+4,1);
% 
% for i=1:NUM_CAMERAS
%     x0((1+7*(i-1)):(7*i)) = camera_calibrations{i}.camera_frame;
% end
% 
% x0(end-3:end-1) = camera_calibrations{1}.marker_position;
% x0(end) = SPHERE_RAD;
% 
% A = zeros(1,7*NUM_CAMERAS+4);
% A(end) = -1;
% options = optimoptions('fmincon','MaxIterations', 4000, 'MaxFunctionEvaluations', 3000);
% x = fmincon(@(x)global_cost(x,trimmed_clouds, robot_poses),x0,[],[],[],[],[],[],@global_const,options);
% 
% for i=1:NUM_CAMERAS
%     camera_calibrations{i}.camera_frame = x((1+7*(i-1)):(7*i))';
%     camera_calibrations{i}.marker_position = x(end-3:end-1);
%     camera_calibrations{i}.sphere_rad = x(end);
% end
% 
% end
% 
% function cost = global_cost(x, trimmed_clouds, robot_poses)
% 
% NUM_CAMERAS = size(trimmed_clouds,1);
% NUM_SAMPLES = size(trimmed_clouds,2);
% 
% camera_poses = zeros(7,NUM_CAMERAS);
% camera_poses(:) = x(1:NUM_CAMERAS*7);
% 
% sphere_center = x((end-3):(end-1));
% sphere_rad = x(end);
% 
% fprintf('%1.2f %1.2f %1.2f %1.2f\n', [sphere_center; sphere_rad]*1000);
% 
% 
% 
% cost = 0;
% for i=1:NUM_CAMERAS
%     frame = camera_poses(:,i)';
%     C = [quat2rotm(frame(4:7)) frame(1:3)'; 0 0 0 1];
%     for j = 1:NUM_SAMPLES
%         if (isempty(trimmed_clouds{i,j}))
%             continue;
%         end
%         frame = robot_poses(j,:);
%         R = [quat2rotm(frame(4:7)) frame(1:3)'; 0 0 0 1];
%         [xc,yc,zc] = nl_sphere(trimmed_clouds{i,j},sphere_rad);
%         diff = R(1:3,1:3) * sphere_center + R(1:3,4) - (C(1:3,1:3) * [xc;yc;zc] + C(1:3,4));
%         cost = cost + diff'*diff;
%     end
% end
% 
% end
% 
% function [c,ceq] = global_const(x)
% 
% c = [];
% 
% NUM_CAMERAS = floor(length(x)/7);
% 
% ceq = zeros(4,1);
% 
% for i=1:NUM_CAMERAS
%     ceq(i) = norm(x((4 + 7*(i-1)):(7*i)));
% end
% 
% 
% end

function [camera_calibrations, trimmed_clouds] = fine_cal(camera_calibrations, all_clouds, robot_poses, SPHERE_RAD, VIEW_ANGLE, PIX_PER_SQRAD)

NUM_CAMERAS = length(camera_calibrations);
NUM_SAMPLES = size(robot_poses,1);

trimmed_clouds = cell(NUM_CAMERAS, NUM_SAMPLES);

% Absolute lower limit of points
MIN_MATCHES = 100;

for i=1:NUM_CAMERAS
    valid_robot_poses = [];
    valid_sphere_centers = [];
    valid_num_matches = [];
    frame = camera_calibrations{i}.camera_frame;
    K = [quat2rotm(frame(4:7)) frame(1:3)'; 0 0 0 1];
    v = camera_calibrations{i}.marker_position';
    for j=1:NUM_SAMPLES
        
        frame = robot_poses(j,:);
        R = [quat2rotm(frame(4:7)) frame(1:3)'; 0 0 0 1];
        
        % Compute where the center of the sphere is according to our
        % current guess of the camera pose
        sphere_center = inv(K)*R*[v;1];
        sphere_center = sphere_center(1:3,1)';
        
        % Only keep points that are within 2 radii of our estimated sphere
        % center
        dist2 = sum(bsxfun(@minus, all_clouds{i}{j}, sphere_center).^2,2);
        keep_idx = (dist2 < (SPHERE_RAD * 2)^2);
        
        % If there are too few points left, we're not going to find a
        % sphere in the remaining points, so move on to the next cloud
        if (sum(keep_idx) < MIN_MATCHES)
            continue;
        end
        
        % Otherwise, fit a sphere to the remaining points
        [xc,yc,zc, num_matches, inliers] = ransac_sphere(all_clouds{i}{j}(keep_idx,:),SPHERE_RAD);
        min_pixels = (2 * atan2(SPHERE_RAD * sin(VIEW_ANGLE), zc - SPHERE_RAD * cos(VIEW_ANGLE)))^2 * PIX_PER_SQRAD;
        if (num_matches > min_pixels)
            % If there were enough matches, save the robot pose and found
            % sphere center
            valid_robot_poses = [valid_robot_poses; robot_poses(j,:)];
            valid_sphere_centers = [valid_sphere_centers; [xc yc zc]];
            valid_num_matches = [valid_num_matches; num_matches];
            % Keep the points that fit the sphere in case we want to do
            % something with them later
            cloud = all_clouds{i}{j}(keep_idx,:);
            trimmed_clouds{i,j} = cloud(inliers,:);
%             [x,y,z] = sphere(30);
%             x = x * SPHERE_RAD + xc;
%             y = y * SPHERE_RAD + yc;
%             z = z * SPHERE_RAD + zc;
%             figure(1);clf;pcshow(cloud(inliers,:), 'r');hold on;
%             pcshow(cloud(setdiff(1:size(cloud,1),inliers),:),'b');
%             lightGrey = 0.8 * [1 1 1];
%             surface(x,y,z,'FaceColor', 'none', 'EdgeColor', lightGrey);
%             drawnow;pause;
        end
    end
    figure(i);clf;histogram(valid_num_matches);drawnow;
    % Now, using all of the found spheres and robot poses, find the pose 
    % of this camera with respect to the robot. Note that here, we set a
    % much higher threshold for the number of poses that have to match,
    % because we're pretty sure we've found a sphere.
    fprintf('num valid poses for camera "%s": %d\n', camera_calibrations{i}.camera_name, size(valid_robot_poses,1));
    [camera_frame, marker_position, average_error, kept_idx] = extrinsic_cal(valid_robot_poses, valid_sphere_centers, 0.9);
    
    % Save all of the relevant info, and move on to the next camera
    camera_calibrations{i}.camera_name = camera_calibrations{i}.camera_name;
    camera_calibrations{i}.camera_frame = camera_frame;
    camera_calibrations{i}.marker_position = marker_position;
    camera_calibrations{i}.average_error = average_error;
    camera_calibrations{i}.robot_poses = valid_robot_poses(kept_idx,:);
    camera_calibrations{i}.camera_positions = valid_sphere_centers(kept_idx,:);
end


end

function [all_clouds, robot_poses, camera_names] = extract_data(folder)
folder_contents = dir(folder);
% First, read in robot poses
robot_poses = [];
NUM_CAMERAS = 0;
camera_names = {};
for i=1:length(folder_contents)
    if (folder_contents(i).name(1) == '.')
        continue;
    end
    if (folder_contents(i).isdir)
        NUM_CAMERAS = NUM_CAMERAS + 1;
        camera_names = [camera_names; folder_contents(i).name];
    else
        robot_poses = load(fullfile(folder,folder_contents(i).name));
    end
end

NUM_SAMPLES = size(robot_poses,1);

all_clouds = cell(NUM_CAMERAS,1);

for i=1:NUM_CAMERAS
    clouds = cell(NUM_SAMPLES,1);
    for j=1:NUM_SAMPLES
        ptcloud = pcread(fullfile(folder, camera_names{i}, sprintf('%d.pcd', j)));
        clouds{j} = ptcloud.Location;
    end
    all_clouds{i} = clouds;
end

end

function [camera_calibrations, all_clouds, robot_poses] = rough_cal(folder, SPHERE_RAD, VIEW_ANGLE, PIX_PER_SQRAD)

% Background setup
% Proportion of images where a voxel needs to be filled to be considered background
BACKGROUND_RATIO = 0.5; 
% Cube in camera's frame that may not be background. Anything outside of 
% this range is considered background, and will be ignored. [xmin xmax;
% ymin ymax; zmin zmax]. Units in meters
% BOUNDING_CUBE = [-0.25 0.25; -0.25 0.25; 0.2 0.8]';
BOUNDING_CUBE = [-1 1; -1 1; 0 1]';
% Side length of a voxel. If this number is too small, everything will be
% considered foreground. If this number is too big, everything will be
% considered background. Units in meters.
BACKGROUND_VOXEL_LENGTH = 0.01;
% The number of voxels in [x,y,z]. Total number of voxels is the product of
% these numbers
BACKGROUND_DIMS = round((BOUNDING_CUBE(2,:) - BOUNDING_CUBE(1,:))/BACKGROUND_VOXEL_LENGTH);

% Absolute lower limit of points
MIN_MATCHES = 50;

folder_contents = dir(folder);
% First, read in robot poses
robot_poses = [];
NUM_CAMERAS = 0;
CAMERA_NAMES = {};
for i=1:length(folder_contents)
    if (folder_contents(i).name(1) == '.')
        continue;
    end
    if (folder_contents(i).isdir)
        NUM_CAMERAS = NUM_CAMERAS + 1;
        CAMERA_NAMES = [CAMERA_NAMES; folder_contents(i).name];
    else
        robot_poses = load(fullfile(folder,folder_contents(i).name));
    end
end

NUM_SAMPLES = size(robot_poses,1);

camera_calibrations = cell(NUM_CAMERAS,1);
all_clouds = cell(NUM_CAMERAS,1);

for i=1:NUM_CAMERAS
    background_count = uint16(zeros(BACKGROUND_DIMS));
    clouds = cell(NUM_SAMPLES,1);
    for j=1:NUM_SAMPLES
        % Currently, no voxels are occupied
        voxel = false(BACKGROUND_DIMS);
        % Get the point cloud
        ptcloud = pcread(fullfile(folder, CAMERA_NAMES{i}, sprintf('%d.pcd', j)));
        clouds{j} = ptcloud.Location;
        
        % Get where each point is between each dimension.
        fracs = bsxfun(@rdivide, bsxfun(@minus, clouds{j}, BOUNDING_CUBE(1,:)), (BOUNDING_CUBE(2,:) - BOUNDING_CUBE(1,:)));
        
        % If any points are outside the bounding cube, remove them
        remove_idx = any(fracs < 0,2) | any(fracs >= 1, 2);
        clouds{j}(remove_idx,:) = [];
        fracs(remove_idx,:) = [];
        
        % Get the voxel indices [x1 y1 z1;x2 y2 z2 ...] 
        idxs = floor(bsxfun(@times, fracs, BACKGROUND_DIMS)) + 1;
        % Now, set the voxels that have a point in them to be true
        voxel(sub2ind(size(voxel),idxs(:,1), idxs(:,2), idxs(:,3))) = true;
        % Add 1 to the background voxels that have a point in them
        background_count(voxel) = background_count(voxel) + 1;
        disp(size(clouds{j}, 1));
    end
    
    valid_robot_poses = [];
    valid_sphere_centers = [];
    % Now, subtract background from each cloud
    for j=1:NUM_SAMPLES
        % Get where each point is between each dimension.
        fracs = bsxfun(@rdivide, bsxfun(@minus, clouds{j}, BOUNDING_CUBE(1,:)), (BOUNDING_CUBE(2,:) - BOUNDING_CUBE(1,:)));
        % Get the voxel indices [x1 y1 z1;x2 y2 z2 ...] 
        idxs = floor(bsxfun(@times, fracs, BACKGROUND_DIMS)) + 1;
        
        % Get the fraction of how often each voxel had a point in it, and
        % if this number is larger than our threshold, it's background
        is_bg = background_count(sub2ind(size(background_count),idxs(:,1), idxs(:,2), idxs(:,3))) > NUM_SAMPLES * BACKGROUND_RATIO;
        
        % Remove all background points
        clouds{j}(is_bg,:) = [];
        
        % If all that's left is less than the minimum number of points
        % we're looking for in a sphere, we're not going to find a sphere
        % in this cloud, so skip to the next cloud.
        if (size(clouds{j},1) < MIN_MATCHES)
            continue;
        end
        
%         % Otherwise, look at the cloud and manually pick depth range
%         pc = pointCloud(clouds{j});
%         figure(1);clf(1);
%         pcshow(pc); zlabel('z'); view(90,0);
%         zmax = input('Max Z:');
%         zmin = input('Min Z:');
%         id_plane = (clouds{j}(:, 3) < zmax) & (clouds{j}(:, 3) > zmin);
%         clouds{j} = clouds{j}(id_plane, :);

        % fit a sphere to the remaining points
        [xc,yc,zc, num_matches] = ransac_sphere(clouds{j},SPHERE_RAD);
        [x,y,z] = sphere(30);
        x = x * SPHERE_RAD + xc;
        y = y * SPHERE_RAD + yc;
        z = z * SPHERE_RAD + zc;
        figure(1);clf;pcshow(clouds{j});hold on;
        lightGrey = 0.8 * [1 1 1];
        surface(x,y,z,'FaceColor', 'none', 'EdgeColor', lightGrey);
        drawnow
        %             drawnow;pause;
        

        % If our sphere fit enough points, remember where the robot and
        % sphere center was
        min_pixels = (2 * atan2(SPHERE_RAD * sin(VIEW_ANGLE), zc - SPHERE_RAD * cos(VIEW_ANGLE)))^2 * PIX_PER_SQRAD;
        if (num_matches > min_pixels)
            valid_robot_poses = [valid_robot_poses; robot_poses(j,:)];
            valid_sphere_centers = [valid_sphere_centers; [xc yc zc]];
            
        end
    end
    % Now, using all of the found spheres and robot poses, find the pose 
    % of this camera with respect to the robot. Note that here, we'll allow
    % for a lower threshold of poses to match, since we still have a lot of
    % points that are probably not the sphere in our images
    fprintf('num valid poses for camera "%s": %d\n', CAMERA_NAMES{i}, size(valid_robot_poses,1));
    [camera_frame, marker_position, average_error] = extrinsic_cal(valid_robot_poses, valid_sphere_centers, 0.5);

    % Save relevant data, and move on to the next camera
    all_clouds{i} = clouds;
    camera_calibrations{i}.camera_name = CAMERA_NAMES{i};
    camera_calibrations{i}.camera_frame = camera_frame;
    camera_calibrations{i}.marker_position = marker_position;
    camera_calibrations{i}.average_error = average_error;
end


end

% Find a sphere of known radius in a point cloud. Uses RANSAC for outlier
% rejection
function [xc,yc,zc, best_matches, best_inliers] = ransac_sphere(p,sphere_rad)
CLOSE_THRESH = sphere_rad / 100;
NUM_ITER = 10000;
N = size(p,1);

best_matches = 0;
best_center = [0;0;0];
best_inliers = [];

for i=1:NUM_ITER
    % First, randomly pick 4 points, and fit a sphere
    [xc,yc,zc] = ls_fit_sphere(p(randperm(N,4),:),sphere_rad);
    
    % Next, compute how far all points are from our sphere, and get a list
    % of points that are close enough to be considered inliers
    dists = sqrt((p(:,1) - xc).^2 + (p(:,2) - yc).^2 + (p(:,3) - zc).^2) - sphere_rad;
    inliers = find(abs(dists) < CLOSE_THRESH);
    
    % If the number of inliers is the best yet, keep this as our best
    % sphere
    num_matches = length(inliers);
    if (num_matches > best_matches)
        best_matches = num_matches;
        best_center = [xc;yc;zc];
        best_inliers = inliers;
    end
end

%fprintf('best_matches: %d\n',best_matches);

% Now, do a nonlinear fit of a sphere with the best points found
[xc,yc,zc] = nl_sphere(p(best_inliers,:), sphere_rad);

end

% Nonlinear fit of a sphere. Uses least-squares fit as an initial guess
function [xc,yc,zc] = nl_sphere(p,r)
p = double(p);
r = double(r);
[xc,yc,zc] = ls_fit_sphere(p,r);
options = optimoptions('fminunc','MaxIterations', 4000, 'MaxFunctionEvaluations', 3000, 'Display', 'off');
x = fminunc(@(x)sphere_cost(x,p,r),[xc;yc;zc],options);
xc = x(1);
yc = x(2);
zc = x(3);
end

% Nonlinear cost function of a sphere
function cost = sphere_cost(x,p,r)
cost = sum(abs(sum(bsxfun(@minus, p, x').^2,2) - r^2));
end

% Least-squares fit of a sphere
function [xc,yc,zc] = ls_fit_sphere(p,r)

N = size(p,1);

A = [2*p ones(N,1)];
b = sum(p.^2,2);

m = pinv(A)*b;

xc = m(1);
yc = m(2);
zc = m(3);

end

