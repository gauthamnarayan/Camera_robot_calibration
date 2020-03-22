% A marker is attached to a robot. The robot is moved to a bunch of
% poses, and a camera finds the 3d position of the marker relative to its
% frame. The position of the marker relative to the robot is unknown, and
% the pose of the camera relative to the robot is unknown. 
% 
% Input: 
% robot poses: [M x 7], where each row is: [x y z q0 qx qy qz] which is the
%       robot's pose (transform from robot base frame to tool frame)
% camera_positions: [M x 3], where each row is [x y z], the 3d position of 
%       the marker in the camera frame.
%
% Output:
% camera_frame: [x y z q0 qx qy qz], The pose of the camera with 
%       respect to the robot's base frame
% marker_position: [x y z], the position of the marker with respect to the
%       robot's tool frame
%
% This function implements the paper:
% Liao Wu and Hongliang Ren "Finding the Kinematic Base Frame of a Robot 
% by Hand-Eye Calibration Using 3D Position Data" in IEEE Transactions on
% Automation Science and Engineering Vol 14, No 1. 2017.
%

function [camera_frame, marker_position, average_error, kept_idx] = extrinsic_cal(robot_poses, camera_positions, ACCEPT_RATIO)

if (nargin == 2)
    ACCEPT_RATIO = 0.9;
end

% Check input
if (size(robot_poses,2) ~= 7 || size(robot_poses,1) ~= size(camera_positions,1) || size(camera_positions,2) ~=3)
    error('robot_poses must be Mx7 and camera_positions must be Mx3');
end

% Compute closed-form solution as a first approximation, with RANSAC
M = size(robot_poses,1);
thresh = norm(max(camera_positions) - min(camera_positions))/36; % reduced threshold from 100

best_err = 100000;
NUM_TRIALS = 3000;
for i=1:NUM_TRIALS
    if (mod(i,100) == 0)
        fprintf('RANSAC trial %d/%d\n',i,NUM_TRIALS);
    end
    idxs = randperm(size(robot_poses,1),10);
    [Rbw, tbw, ph] = closed_form_approximation(robot_poses(idxs,:), camera_positions(idxs,:));
    [Rbw, tbw, ph] = iterative_solution(robot_poses(idxs,:), camera_positions(idxs,:), Rbw, tbw, ph);
    
    Hbw = [Rbw tbw; 0 0 0 1];
    marker_positions = zeros(4,M);
    for j=1:M
        Hbh = [quat2rot(robot_poses(j,4:7)) robot_poses(j,1:3)'; 0 0 0 1];
        marker_positions(:,j) = inv(Hbh) * Hbw * [camera_positions(j,:)';1];
    end
    errors = sqrt(sum(bsxfun(@minus, marker_positions(1:3,:), ph).^2));
    keep_idx = (errors < thresh);
    
    if(sum(keep_idx) < ACCEPT_RATIO * M)
        continue;
    end
    
    [Rbw, tbw, ph] = closed_form_approximation(robot_poses(keep_idx,:), camera_positions(keep_idx,:));
    [Rbw, tbw, ph] = iterative_solution(robot_poses(keep_idx,:), camera_positions(keep_idx,:), Rbw, tbw, ph);
    idxs = find(keep_idx);
    Hbw = [Rbw tbw; 0 0 0 1];
    marker_positions = zeros(4,length(idxs));
    for j= 1:length(idxs)
      Hbh = [quat2rot(robot_poses(idxs(j),4:7)) robot_poses(idxs(j),1:3)'; 0 0 0 1];
      marker_positions(:,j) = inv(Hbh) * Hbw * [camera_positions(idxs(j),:)';1];
    end
    errors = sqrt(sum(bsxfun(@minus, marker_positions(1:3,:), ph).^2));
    score = mean(errors);
    if (score < best_err)
        best_err = score;
        best_idx = keep_idx;
        best_Rbw = Rbw;
        best_tbw = tbw;
        best_ph = ph;
    end
end

fprintf('Poses kept: %d/%d\n',sum(best_idx),M);

% Form our output vectors
camera_frame = [best_tbw', rot2quat(best_Rbw)'];
marker_position = best_ph';
average_error = best_err;
kept_idx = best_idx;

end



function [Rbw, tbw, ph] = closed_form_approximation(robot_poses, camera_positions)
% Get number of poses
M = size(robot_poses,1);

% First, preallocate matrix sizes for the equation K * x = b
K = zeros(3*M,15);
b = zeros(3*M,1);
for i=1:M
    Rbh = quat2rot(robot_poses(i,4:7));
    
    K((3*(i-1)+1):(3*i),:) = [kron(camera_positions(i,:),eye(3)), eye(3), -1*Rbh];
    b((3*(i-1)+1):(3*i)) = robot_poses(i,1:3);
end

% Now, compute the least squares solution. x contains the following: the
% rotation matrix and translation from the robot base frame to the camera
% frame, and the position of the marker in the hand frame, which are our
% unknowns
x = pinv(K)*b;

% Extract the rotation matrix approximation, and orthogonalize it
Rbwt = [x(1:3,1), x(4:6,1), x(7:9,1)];
[u,s,v] = svd(Rbwt);

Rbw = sign(det(s))*u*v';

% Extract the translational component and marker position
tbw = x(10:12,1);
ph = x(13:15,1);
end


function [Rbw, tbw, ph] = iterative_solution(robot_poses, camera_positions, Rbw, tbw, ph)
% Get number of poses
M = size(robot_poses,1);

% Preallocate jacobian and cost function
J=zeros(3*M, 9);
f0 = zeros(3*M,1); 

count = 0;
while count < 1000
    % Convert matrix into skew symmetric matrix
    w = rot2skew(Rbw);
    % Compute derivative of Rbw
    if (w == 0)
        Qbw = Rbw';
    else
        Qbw = Rbw' * (eye(3) + (1-cos(norm(w)))/(norm(w))^2*skew_sym(w) + (norm(w)-sin(norm(w)))/(norm(w)^3)*skew_sym(w)^2);
    end
    % Create Jacobian and cost function matrices
    for i=1:M
        Rbh = quat2rot(robot_poses(i,4:7));
        pw = camera_positions(i,:)';
        tbh = robot_poses(i,1:3)';
        
        J((3*(i-1)+1):(3*i),:) = [-Rbw*skew_sym(pw)*Qbw, eye(3), -Rbh];
        f0((3*(i-1)+1):(3*i)) = Rbw * pw + tbw - Rbh * ph - tbh;
    end
    % Compute our step
    y = -pinv(J)*f0;
    
    % Update our matrices with our step
    Rbw = skew2rot(w + y(1:3,1));
    tbw = tbw + y(4:6,1);
    ph = ph + y(7:9,1);
    
    % If step size was small enough, we're done!
    if (norm(y) < 0.00001)
        break;
    end
    
    count = count + 1;
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

