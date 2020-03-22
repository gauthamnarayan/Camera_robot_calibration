import rospy
#import sys
#import pypcd
#from sensor_msgs.msg import PointCloud2
#from sensor_msgs.msg import Image
#import pprint
#import cv2
#from cv_bridge import CvBridge
import tf
import time

def main():
	
	listener = tf.TransformListener()
 	rospy.init_node("test_listener")
	
	while not rospy.is_shutdown():
		try:
                        print("trying")
			trans,rot = listener.lookupTransform('/reference/head', '/reference/right_l6', rospy.Time(0))
			print(trans)
			print(rot)
		except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException) as e:
			print(e)
			continue


if __name__ == '__main__':
	main()




