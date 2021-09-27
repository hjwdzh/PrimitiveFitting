#include <opencv2/opencv.hpp>
#include "PlaneFitter.h"
#include <glm/glm.hpp>
#include <vector>
#include <map>
#include "Plane.h"

std::pair<cv::Mat, std::vector<PlaneHelper> > ExtractPlane(cv::Mat depth_input, double fx, double fy, double cx, double cy, double threshold)
{
	cv::Mat depth;
	if (depth_input.type() == CV_16UC1) {
		depth_input.convertTo(depth, CV_32F, 1e-3);
	}
	else {
		depth = depth_input.clone();
	}

	PlaneFitter pf;
	PointCloud cloud(depth, fx, fy, cx, cy);
	cv::Mat seg(depth.rows, depth.cols, CV_8UC3);
	pf.Process(&cloud);

	std::vector<std::vector<glm::vec3> > points_list;
	
	std::map<int, int> type_list;
	std::vector<float> scores;
	for (int i = 0; i < pf.tags.rows; ++i) {
		for (int j = 0; j < pf.tags.cols; ++j) {
			int tag = pf.tags.at<int>(i, j);
			if (tag < 0)
				continue;
			auto it = type_list.find(tag);
			int idx = type_list.size();
			if (it == type_list.end()) {
				type_list.insert(std::make_pair(tag, type_list.size()));
				points_list.push_back(std::vector<glm::vec3>());
			}
			else {
				idx = it->second;
			}
			double x, y, z;
			if (cloud.get(i, j, x, y, z)) {
				points_list[idx].push_back(glm::vec3(x, y, z));
			}
		}
	}
	scores.resize(points_list.size());
	std::vector<PlaneHelper> planeHelpers(points_list.size());
	std::vector<int> compact_indices(points_list.size());
	int offset = 0;
	for (int i = 0; i < points_list.size(); ++i) {
		glm::vec3 center(0);
		for (int j = 0; j < points_list[i].size(); ++j)
			center += points_list[i][j];
		center /= points_list[i].size();
		cv::Mat A(3, 3, CV_64F);
		memset(A.data, 0, sizeof(double) * 9);
		memset(planeHelpers[i].disMat, 0, sizeof(double) * 16);
		for (int j = 0; j < points_list[i].size(); ++j) {
			glm::vec3 diff = 1e3f * (points_list[i][j] - center);
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					A.at<double>(k, l) += diff[k] * diff[l];
				}
			}
			double w[4];
			w[0] = points_list[i][j].x;
			w[1] = points_list[i][j].y;
			w[2] = points_list[i][j].z;
			w[3] = 1;
			
			for (int i1 = 0; i1 < 4; ++i1) {
				for (int j1 = 0; j1 < 4; ++j1) {
					planeHelpers[i].disMat[i1 * 4 + j1] += w[i1] * w[j1];
				}
			}
		}


		A /= points_list[i].size();
		cv::SVD thissvd(A, cv::SVD::FULL_UV);

		cv::Mat S = thissvd.w;

		scores[i] = S.at<double>(2, 0) / S.at<double>(1, 0);
		if (scores[i] > threshold)
			compact_indices[i] = -1;
		else
			compact_indices[i] = offset++;
		double len = sqrt(thissvd.vt.at<double>(2, 0)*thissvd.vt.at<double>(2, 0) + thissvd.vt.at<double>(2, 1)*thissvd.vt.at<double>(2, 1)
			+ thissvd.vt.at<double>(2, 2)*thissvd.vt.at<double>(2, 2));
		planeHelpers[i].params[0] = thissvd.vt.at<double>(2, 0) / len;
		planeHelpers[i].params[1] = thissvd.vt.at<double>(2, 1) / len;
		planeHelpers[i].params[2] = thissvd.vt.at<double>(2, 2) / len;
		planeHelpers[i].params[3] = (-center.x * thissvd.vt.at<double>(2, 0) - center.y * thissvd.vt.at<double>(2, 1) - center.z * thissvd.vt.at<double>(2, 2))/len;
	}
	for (int i = 0; i < compact_indices.size(); ++i) {
		if (compact_indices[i] >= 0)
			planeHelpers[compact_indices[i]] = planeHelpers[i];
	}
	planeHelpers.resize(offset);
	cv::Mat tags = cv::Mat(pf.tags.size(), CV_8U);
	memset(tags.data, 0, sizeof(unsigned char) * tags.cols * tags.rows);
	for (int i = 0; i < pf.tags.rows; ++i) {
		for (int j = 0; j < pf.tags.cols; ++j) {
			int tag = pf.tags.at<int>(i, j);
			if (tag < 0) {
				tags.at<unsigned char>(i, j) = -1;
				continue;
			}
			auto it = type_list.find(tag);
			int idx = it->second;
			if (compact_indices[idx] >= 0) {
				tags.at<unsigned char>(i, j) = compact_indices[idx];
			}
			else {
				tags.at<unsigned char>(i, j) = -1;
			}
		}
	}

	return std::make_pair(tags, planeHelpers);
}