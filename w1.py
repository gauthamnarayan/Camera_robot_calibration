import rospy
import sys
import pypcd
from sensor_msgs.msg import PointCloud2
from sensor_msgs.msg import Image
import pprint
import cv2
from cv_bridge import CvBridge
import tf
import time

import moveit_commander
import moveit_msgs.msg
import geometry_msgs.msg
from math import pi
from std_msgs.msg import String
from geometry_msgs.msg import Pose, Point, Quaternion
from moveit_commander.conversions import pose_to_list
import intera_interface
import quat_math

global saved_PCD
global saved_Image
global fileID
saved_PCD = False
saved_Image = False
saved_Pose = False
captureTime = None

def callbackWritePCD(data):
    global saved_PCD, saved_Image, saved_Pose, fileID
    if(saved_Image==True and saved_PCD==True and saved_Pose==True): 
        rospy.signal_shutdown("")

    pc = pypcd.PointCloud.from_msg(data)
    if(saved_PCD == False):
        #pc.save_pcd("pointcloud.pcd", compression='binary')
        pc.save_pcd("./raw_data/" + fileID + ".pcd", compression='binary')
        saved_PCD = True
        print("PCD Saved")
        sys.exit()

def callbackWriteImage(data):
    global saved_PCD, saved_Image, saved_Pose, fileID
    if(saved_Image==True and saved_PCD==True and saved_Pose==True): 
        rospy.signal_shutdown("")

    bridge = CvBridge()
    cv_image = bridge.imgmsg_to_cv2(data, desired_encoding="passthrough")
    cv_image = cv2.cvtColor(cv_image, cv2.COLOR_BGR2RGB)
    if(saved_Image == False):
        global captureTime
        captureTime = data.header.stamp
        cv2.imwrite("./raw_data/" + fileID + '.jpg', cv_image)
        print("Image Saved")
        saved_Image = True
        sys.exit()

def callbackWritePose(data):
    global captureTime

    trans,rot = listener.lookupTransform('base', 'right_hand')

#    print(captureTime.secs)
#    global saved_PCD, saved_Image, saved_Pose, fileID
#    if(saved_Image==True and saved_PCD==True and saved_Pose==True): 
#        rospy.signal_shutdown("")
#
#    if(saved_Pose == False):
#        endpointPose = limb.endpoint_pose()
#	position = endpointPose['position']
#	cartX = position.x
#	cartY = position.y
#	cartZ = position.z
#	orientation = endpointPose['orientation']
#	quatX = orientation.x
#	quatY = orientation.y
#	quatZ = orientation.z
#	quatW = orientation.w
#
#	print("%s,%s,%s,%s,%s,%s,%s\n" %(cartX,cartY,cartZ,quatW,quatX,quatY,quatZ))
#
#        #posesFile = open("./raw_data/robot_poses.txt", 'a')
#	#posesFile.write("%s,%s,%s,%s,%s,%s,%s\n" %(cartX,cartY,cartZ,quatW,quatX,quatY,quatZ))
#	#posesFile.close() 
#        #print(endpointPose)
#        saved_Pose = True


#fileID = sys.argv[1]
#
rospy.init_node('capture_pose_image_pcd')
#
#moveit_commander.roscpp_initialize(sys.argv)
#scene = moveit_commander.PlanningSceneInterface()
#group = moveit_commander.MoveGroupCommander("right_arm")
#limb = intera_interface.Limb("right")

# These are the subscriber topics for the Azure Kinect
rospy.Subscriber("/points2", PointCloud2, callbackWritePCD)
rospy.Subscriber("/rgb/image_raw", Image, callbackWriteImage)
rospy.Subscriber("/rgb/image_raw", Image, callbackWritePose)

rospy.spin()




