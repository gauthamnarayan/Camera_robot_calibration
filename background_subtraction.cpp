#include <iostream>
#include <stdlib.h>

#include <opencv2/highgui.hpp>
#include <opencv2/aruco.hpp>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/passthrough.h>

using namespace std;
using namespace cv;
using namespace pcl;

int main(int argc, char **argv)
{
	//Input 1: drawImage(yes/no)
    //Input 2: folderPath
    //Input 3: inputNumber
    //Input 4: BBLimit
    //Input 5: AR Tag side length in meters

	string drawImage = argv[1];
    string folderPath = argv[2];
    string inputNumber = argv[3];
    string inputImage = folderPath + "/raw_data/" + inputNumber + ".jpg";
    string inputCloud = folderPath + "/raw_data/" + inputNumber + ".pcd";
    string outputCloud = folderPath + "/cropped_pcd/" + inputNumber + ".pcd";
    double BBLimit = (double) atof(argv[4]);
    //double tagLength = (double) atof(argv[5]);

    // Read image using an OpenCV function.
    Mat image = cv::imread(inputImage, IMREAD_COLOR);
    cout<<"Image width: "<<image.size().width<<endl;
    cout<<"Image height!: "<<image.size().height<<endl;
    
    // Set April tag dictionary and create detector parameters.
    Ptr<aruco::DetectorParameters> detectorParams = aruco::DetectorParameters::create();
    Ptr<aruco::Dictionary> dictionary = aruco::getPredefinedDictionary(cv::aruco::DICT_6X6_250);

    // Detect April Tag corners.
    vector< int > ids;
    vector< vector< Point2f > > corners, rejected;
    vector< Vec3d > rvecs, tvecs;
    Vec3d vec;
    aruco::detectMarkers(image, dictionary, corners, ids, detectorParams, rejected);
    
	aruco::drawDetectedMarkers(image, corners, ids);

    // Estimate pose of April Tag using camera matrix & Distortion coeff's - (RosCameraInfo of AzureKinect node).
    double distCoeffArray[8] = {0.6862431764602661, -2.8917713165283203, 0.0004992862232029438, -4.462565993890166e-05, 1.6113708019256592, 0.5638872385025024, -2.7169768810272217, 1.540696382522583};
	// Camera matrix for Intel RealSense
    //double camArray[3][3] = {{616.0628662109375, 0.0, 325.2831726074219},
    //                         {0.0, 616.201416015625, 234.31954956054688},
    //                         {0.0, 0.0, 1.0}};
	// Camera matrix for AzureKinect
    double camArray[3][3] = {{611.1693115234375, 0.0, 638.8131713867188},
                             {0.0, 611.1657104492188, 367.74871826171875},
                             {0.0, 0.0, 1.0}};
    Mat distCoeffMat(1, 8, CV_64F,distCoeffArray);
    Mat cameraMat(3,3,CV_64F,camArray); //Keep double precision to avoid pointer errors.

    // 0.04 is the side length of the QRTag in Meters.
    aruco::estimatePoseSingleMarkers(corners, 0.04, cameraMat, distCoeffMat, rvecs, tvecs);

    for(int i=0; i<ids.size(); i++)
    {
        aruco::drawAxis(image, cameraMat, distCoeffMat, rvecs[i], tvecs[i], 0.04);
    }

    if(drawImage=="yes")
    {
        imshow("", image);
        waitKey(0);
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_filtered (new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_filtered_nonan (new pcl::PointCloud<pcl::PointXYZRGB>);

    if (pcl::io::loadPCDFile<pcl::PointXYZRGB> (inputCloud, *cloud) == -1) //* load the file
    {   
        PCL_ERROR ("Couldn't read input pcd file \n");
        return (-1);
    }

    pcl::PointXYZRGB p = (*cloud)[corners[0][0].y * cloud->width + corners[0][0].x];
    vec[0] = p.x;
    vec[1] = p.y;
    vec[2] = p.z;
    
   double BBSize = BBLimit/2; // In Meters 
   double xLimits[] = {vec[0]-BBSize, vec[0]+BBSize};
   double yLimits[] = {vec[1]-BBSize, vec[1]+BBSize};
   double zLimits[] = {vec[2]-BBSize, vec[2]+BBSize};
  
   // Use below command based on PCL version.	
   //pcl::io::copyPointCloud(*cloud, *cloud_filtered);
    pcl::copyPointCloud(*cloud, *cloud_filtered);
   
    pcl::PassThrough<pcl::PointXYZRGB> passX;
    passX.setInputCloud (cloud_filtered);
    passX.setFilterFieldName ("x");
    passX.setFilterLimits (xLimits[0], xLimits[1]);
    passX.filter (*cloud_filtered);

    pcl::PassThrough<pcl::PointXYZRGB> passY;
    passX.setInputCloud (cloud_filtered);
    passX.setFilterFieldName ("y");
    passX.setFilterLimits (yLimits[0], yLimits[1]);
    passX.filter (*cloud_filtered);

    pcl::PassThrough<pcl::PointXYZRGB> passZ;
    passX.setInputCloud (cloud_filtered);
    passX.setFilterFieldName ("z");
    passX.setFilterLimits (zLimits[0], zLimits[1]);
    passX.filter (*cloud_filtered);

    std::vector<int> index;
    pcl::removeNaNFromPointCloud(*cloud_filtered, *cloud_filtered_nonan, index);
    
	if(pcl::io::savePCDFileASCII (outputCloud, *cloud_filtered_nonan) == 0)
	{
		cout<<"Output pointcloud: "<<outputCloud<<endl;
	}
	else{
		cout<<"File not written"<<endl;
	}
        
    return 0;
}

