# Computing camera extrinsics:

The idea is to provide a known 3d pose of an object in the robot's reference frame and the position of the same object in the camera's reference frame. We choose a cue ball to be our object as it is possible to calculate the center of the sphere from various view points and ball orientations using RANSAC.

Our current procedure for camera-robot calibration:

1. Collect 3d point-clouds of the cue ball attached to the robot arm and also the pose of the robot link attached to the ball(right_l6 for Sawyer). 

2. Remove points that we know are far away from the ball. We want to keep nearby points to improve the precision of the RANSAC sphere fit. We use an April tag to detect the 3d location of the cue ball in the pointcloud and remove all points beyond a threshold distance from the ball.

3. Provide this data to a matlab script(Courtesy-Robbie Paolini) to compute the RANSAC sphere fit and also optimize the camera extrinsic pose over object poses provided.

Usually ~10 poses should be sufficient to estimate the camera extrinsics. The script also produces a visualisation of the sphere fit, this is usefull to debug isssues.

Run main.m to with the following folder structure to obtain the camera extrinsic parameters. 
.  
├── camera_1  
│   ├── 10.pcd  
│   ├── 1.pcd  
│   ├── 2.pcd  
│   ├── 3.pcd  
│   ├── 4.pcd  
│   ├── 5.pcd  
│   ├── 6.pcd  
│   ├── 7.pcd  
│   ├── 8.pcd  
│   └── 9.pcd  
└── robot_poses.txt  

robot_poses.txt is a N x 7 matrix of [x y z q0 qx qy qz]
  
# Background points removal:

Compile the background subtraction code from [here](https://github.com/gauthamnarayan/Camera_robot_calibration/blob/master/background_subtraction.cpp):   
Dependencies: PCL and OpenCV C++ Libraries

Run the code as   
./background_subtraction  drawImage(yes/no) folderPath  inputFileNumber(1/2/3...10) BBLimit Size_of_AprilTag(meters)  


