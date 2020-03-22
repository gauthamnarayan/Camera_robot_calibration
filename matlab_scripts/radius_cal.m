function [camera_frames, marker_position, new_camera_positions, average_errors] = radius_cal(robot_poses, camera_positions)

options = optimoptions('fminunc','MaxIterations', 4000, 'MaxFunctionEvaluations', 3000);
dr = fminunc(@(dr)radius_cost(robot_poses, camera_positions,dr),0,options);

dr

new_camera_positions = correct_camera_positions(camera_positions, dr);

[camera_frames, marker_position, average_errors] = multi_extrinsic_cal(robot_poses, new_camera_positions);

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
    corr = (dists + dr).*(dists+dr > 0) + dists .* (dists+dr <= 0);
    new_camera_positions{i} = bsxfun(@times, bsxfun(@rdivide, camera_positions{i}, dists), corr); 
end
end