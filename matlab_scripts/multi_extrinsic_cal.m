% A marker is attached to a robot. The robot is moved to a bunch of
% poses, and N cameras find the 3d position of the marker relative to its
% frame. The position of the marker relative to the robot is unknown, and
% the pose of each camera relative to the robot is unknown. This function
% takes in the robot positions and found marker positions in each camera
% frame, and returns the frames of each of the cameras with respect to the
% robot base frame.
% 
% Input: 
% robot poses: A Nx1 cell array of [M_n x 7], where each row is: 
%       [x y z q0 qx qy qz] which is the robot's pose (transform from 
%       robot base frame to tool frame) for each of the N cameras
% camera_positions: A Nx1 cell array of [M_n x 3], where each row is 
%       [x_1 y_1 z_1; x_2 y_2 z_2; x_N, y_N, z_N], the 3d position of 
%       the marker in each of the camera frames.
%
% Output:
% camera_frames: a Nx7 matrix, where each row is [x y z q0 qx qy qz], 
%       the pose of each camera with respect to the robot's base frame
% marker_position: [x y z], the position of the marker with respect to the
%       robot's tool frame
% average_errors: a Nx1 array with the average error of the 3d marker
%       position for each camera over all of the images.
%
% This function extends the paper:
% Liao Wu and Hongliang Ren "Finding the Kinematic Base Frame of a Robot 
% by Hand-Eye Calibration Using 3D Position Data" in IEEE Transactions on
% Automation Science and Engineering Vol 14, No 1. 2017.
%
% By applying it to multiple cameras looking at a single target.
%

function [camera_frames, marker_position, average_errors] = multi_extrinsic_cal(robot_poses, camera_positions)

% Check input
if (~iscell(robot_poses) || ~iscell(camera_positions) || length(robot_poses) ~= length(camera_positions))
    error('robot_poses and camera_positions must both be cell arrays of the same length');
else
    for i=1:length(robot_poses)
        rp = robot_poses{i};
        cp = camera_positions{i};
        if (size(rp,2) ~= 7 || size(rp,1) ~= size(cp,1) || size(cp,2) ~= 3)
            error('for camera %d robot_poses must be Mx7 and camera_positions must be Mx3', i);
        end
    end
end

N = length(robot_poses);

[Rbw, tbw, ph] = closed_form_approximation(robot_poses, camera_positions);

% Compute average error
% average_errors = zeros(N,1);
% for i=1:N
%     robot_pose = robot_poses{i};
%     camera_position = camera_positions{i};
%     err_sum = 0;
%     for j=1:size(robot_pose,1)
%         err_sum = err_sum + norm(Rbw{i} * camera_position(j,:)' + tbw{i} - (quat2rot(robot_pose(j,4:7)) * ph + robot_pose(j,1:3)'));
%     end
%     average_errors(i) = err_sum / size(robot_pose,1);
% end
% 
% fprintf('ph = [%1.4f, %1.4f, %1.4f]\n', ph);
% for i=1:N
%     fprintf('%1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f (%1.2f mm)\n', tbw{i}, rot2quat(Rbw{i}), average_errors(i) * 1000);
% end
% fprintf('\n');




[Rbw, tbw, ph] = iterative_solution(robot_poses, camera_positions, Rbw, tbw, ph);

% Compute average error
average_errors = zeros(N,1);
for i=1:N
    robot_pose = robot_poses{i};
    camera_position = camera_positions{i};
    err_sum = 0;
%     errors = zeros(size(robot_pose,1),1);
    for j=1:size(robot_pose,1)
%         errors(j) = norm(Rbw{i} * camera_position(j,:)' + tbw{i} - (quat2rot(robot_pose(j,4:7)) * ph + robot_pose(j,1:3)'));
        err_sum = err_sum + norm(Rbw{i} * camera_position(j,:)' + tbw{i} - (quat2rot(robot_pose(j,4:7)) * ph + robot_pose(j,1:3)'));
    end
    average_errors(i) = err_sum / size(robot_pose,1);
%     figure(i);clf;hist(errors,10);
end

% fprintf('ph = [%1.4f, %1.4f, %1.4f]\n', ph);
% for i=1:N
%     fprintf('%1.4f, %1.4f, %1.4f, %1.5f, %1.5f, %1.5f, %1.5f (%1.2f mm)\n', tbw{i}, rot2quat(Rbw{i}), average_errors(i) * 1000);
% end



camera_frames = cell(N,1);
% Form our output vectors
for i=1:N
    camera_frames{i} = [tbw{i}', rot2quat(Rbw{i})'];
end
marker_position = ph';

end



function [Rbw, tbw, ph] = closed_form_approximation(robot_poses, camera_positions)
% Get number of cameras
N = length(robot_poses);

% We will form large matrices Kt * x = bt
% Kt is a large stacked matrix
Kt = [];
bt = [];
for j=1:N
    M = size(robot_poses{j},1);
    K = zeros(3*M,12*N+3);
    b = zeros(3*M,1);
    for i=1:M
        Rbh = quat2rot(robot_poses{j}(i,4:7));

        K((3*(i-1)+1):(3*i),:) = [zeros(3,12*(j-1)), kron(camera_positions{j}(i,:),eye(3)), eye(3), zeros(3,12*(N-j)), -1*Rbh];
        b((3*(i-1)+1):(3*i)) = robot_poses{j}(i,1:3);
    end
    Kt = [Kt;K];
    bt = [bt;b];
end



% Now, compute the least squares solution. x contains the following: the
% rotation matrix and translation from the robot base frame to the camera
% frame, for each camera frame, and the position of the marker in the 
% hand frame, which are our unknowns
x = pinv(Kt)*bt;

Rbw = cell(N,1);
tbw = cell(N,1);
for j=1:N
    % Extract the rotation matrix approximation, and orthogonalize it
    Rbwt = [x((12*(j-1)+1):(12*(j-1)+3),1), x((12*(j-1)+4):(12*(j-1)+6),1), x((12*(j-1)+7):(12*(j-1)+9),1)];
    [u,s,v] = svd(Rbwt);
    Rbw{j} = sign(det(s))*u*v';

    % Extract the translational component
    tbw{j} = x((12*(j-1)+10):(12*(j-1)+12),1);
end

% Extract marker position
ph = x(end-2:end,1);

end


function [Rbw, tbw, ph] = iterative_solution(robot_poses, camera_positions, Rbw, tbw, ph)
% Get number of cameras
N = length(robot_poses);

wc = cell(N,1);
MAX_STEPS = 100;
for steps = 1:MAX_STEPS
    Jt = [];
    f0t = [];
    for j=1:N
        % Get number of poses
        M = size(robot_poses{j},1);

        % Preallocate jacobian and cost function
        J=zeros(3*M, 6*N+3);
        f0 = zeros(3*M,1); 
    
        % Convert matrix into skew symmetric matrix
        w = rot2skew(Rbw{j});
        wc{j} = w;
        % Compute derivative of Rbw
        if (w == 0)
            Qbw = Rbw{j}';
        else
            Qbw = Rbw{j}' * (eye(3) + (1-cos(norm(w)))/(norm(w))^2*skew_sym(w) + (norm(w)-sin(norm(w)))/(norm(w)^3)*skew_sym(w)^2);
        end
        % Create Jacobian and cost function matrices
        for i=1:M
            Rbh = quat2rot(robot_poses{j}(i,4:7));
            pw = camera_positions{j}(i,:)';
            tbh = robot_poses{j}(i,1:3)';

            J((3*(i-1)+1):(3*i),:) = [zeros(3,6*(j-1)), -Rbw{j}*skew_sym(pw)*Qbw, eye(3), zeros(3,6*(N-j)), -Rbh];
            f0((3*(i-1)+1):(3*i)) = Rbw{j} * pw + tbw{j} - Rbh * ph - tbh;
        end
        
        Jt = [Jt;J];
        f0t = [f0t;f0];
    end
    % Compute our step
    y = -pinv(Jt)*f0t;
    
    % Update our matrices with our step
    for j=1:N
        Rbw{j} = skew2rot(wc{j} + y((6*(j-1)+1):(6*(j-1)+3),1));
        tbw{j} = tbw{j} + y((6*(j-1)+4):(6*(j-1)+6),1);
    end
    ph = ph + y(end-2:end,1);
    
    % If step size was small enough, we're done!
    if (norm(y) < 0.000000001)
        break;
    end
end

if (norm(y) >= 0.000000001)
    warning('maximum iterations reached. Unknown behavior');
end

end

function w = rot2skew(R)
theta = acos((trace(R)-1)/2);
if (theta == 0)
    w = [0;0;0];
else
    w = theta/(2*sin(theta))*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
end
end

function R = skew2rot(w)
R = eye(3) + sin(norm(w))/norm(w)*skew_sym(w) + (1-cos(norm(w)))/(norm(w)^2)*skew_sym(w)^2;
end

function w_hat = skew_sym(w)
w_hat = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

function R = quat2rot(q)

q = q/norm(q);

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

R =[q0^2+q1^2-q2^2-q3^2,    2*(q1*q2-q0*q3),        2*(q1*q3+q0*q2);
    2*(q1*q2+q0*q3),        q0^2-q1^2+q2^2-q3^2,    2*(q2*q3-q0*q1);
    2*(q1*q3-q0*q2),        2*(q2*q3+q0*q1),        q0^2-q1^2-q2^2+q3^2];
end

% Convert rotation matrix into quaternion
function q=rot2quat(R)

q2 = zeros(4,1);

q2(1) = 1/4*(1+R(1,1)+R(2,2)+R(3,3));
q2(2) = 1/4*(1+R(1,1)-R(2,2)-R(3,3));
q2(3) = 1/4*(1-R(1,1)+R(2,2)-R(3,3));
q2(4) = 1/4*(1-R(1,1)-R(2,2)+R(3,3));

[~,i] = max(q2);

switch i
    case 1
        q0 = sqrt(q2(1));
        q1 = 1/(4*q0)*(R(3,2)-R(2,3));
        q2 = 1/(4*q0)*(R(1,3)-R(3,1));
        q3 = 1/(4*q0)*(R(2,1)-R(1,2));
    case 2
        q1 = sqrt(q2(2));
        q0 = 1/(4*q1)*(R(3,2)-R(2,3));
        q2 = 1/(4*q1)*(R(1,2)+R(2,1)); 
        q3 = 1/(4*q1)*(R(1,3)+R(3,1)); 
    case 3
        q2 = sqrt(q2(3));
        q0 = 1/(4*q2)*(R(1,3)-R(3,1));
        q1 = 1/(4*q2)*(R(1,2)+R(2,1));
        q3 = 1/(4*q2)*(R(2,3)+R(3,2));
    case 4
        q3 = sqrt(q2(4));
        q0 = 1/(4*q3)*(R(2,1)-R(1,2));
        q1 = 1/(4*q3)*(R(1,3)+R(3,1));
        q2 = 1/(4*q3)*(R(2,3)+R(3,2));
    otherwise
        error('something went wrong');
end

q = [q0;q1;q2;q3];

end