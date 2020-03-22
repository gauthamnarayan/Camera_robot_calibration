function error_matrix = rs_check_calibration(camera_frames, marker_position, robot_poses, camera_positions)

NUM_CAMERAS = length(camera_frames);
error_matrix = zeros(NUM_CAMERAS);
for i=1:NUM_CAMERAS
    c = camera_frames{i};
    R = quat2rotm(c(4:7));
    t = c(1:3)';
    for j=(i+1):NUM_CAMERAS
        
        % Get the poses that were common to both cameras
        joint_robot_poses = [];
        camera_positions1 = [];
        camera_positions2 = [];
        for k=1:size(robot_poses{i},1)
            for m=1:size(robot_poses{j},1)
                if all(robot_poses{i}(k,:) == robot_poses{j}(m,:))
                    joint_robot_poses = [joint_robot_poses;robot_poses{i}(k,:)];
                    camera_positions1 = [camera_positions1; camera_positions{i}(k,:)];
                    camera_positions2 = [camera_positions2; camera_positions{j}(m,:)];
                end
            end
        end
        
        % Now compute the error
        c = camera_frames{j};
        R2 = quat2rotm(c(4:7));
        t2 = c(1:3)';
        
        error_matrix(i,j) = mean(sqrt(sum((bsxfun(@plus, R * camera_positions1', t) - bsxfun(@plus, R2 * camera_positions2', t2)).^2)));
    end
end

error_matrix = error_matrix + error_matrix';




end

