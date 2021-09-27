#ifndef PLANE_H_
#define PLANE_H_

#ifndef PLANE_LIB
#define EXPORT_H
#else
#define EXPORT_H __declspec(dllimport)
#endif

#include <opencv2/opencv.hpp>

class PlaneHelper
{
public:
	double params[4];
	double disMat[16];
};

EXPORT_H std::pair<cv::Mat, std::vector<PlaneHelper> > ExtractPlane(cv::Mat depth_input, double fx, double fy, double cx, double cy, double threshold = 1e30);

#endif