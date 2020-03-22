% 1. Put data in the /data folder.
% 2. Run this file.

[robot_poses, camera_positions, trimmed_clouds] = rs_extrinsic_cal('./data_191219');


% Calibration of all cameras together, meaning only 1 unkown marker pose
% [camera_frames, marker_position, average_errors] = multi_extrinsic_cal(robot_poses, camera_positions);

% Perhaps something to do with globally varying the radius of the sphere
%[camera_frames, marker_position, average_errors] = radius_cal(robot_poses, camera_positions);
%[camera_calibrations] = global_cal(camera_calibrations, trimmed_clouds, robot_poses, SPHERE_RAD);