#ifndef PLAnumE_FITTER_H_
#define PLAnumE_FITTER_H_

#include <vector>
#include <set>
#include <queue>
#include <map>
#define _USE_MATH_DEFInumES
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "PlaneSegment.h"

class PlaneFitter {
public:
	typedef std::priority_queue<PlaneSegment::shared_ptr,
		std::vector<PlaneSegment::shared_ptr>,
		MSEHelper> QueueMSE;

	PlaneFitter();
	~PlaneFitter();
	void clear();
	int Neighbors(const int i, const int j, const int H, const int W, int neighbors[4]);
	int BlockID(const int cx, const int cy) const;
	void InitGraph(QueueMSE& queueMSE);
	void MergeGraph(QueueMSE& queueMSE);
	void Refine();
	void AnalyzeBlockTags(std::vector<bool>& valid);
	void FloodFill();
	void Process(const PointCloud* _points);

	const PointCloud *points;
	int width, height;

	int max_group_iter;
	int min_block_size;
	int windowWidth;
	int windowHeight;

	Parameters params;

	shared_ptr<DisjointSet> disjoint_set;
	std::vector<PlaneSegment::shared_ptr> planes;
	cv::Mat tags;

	std::map<int, int> root_block;
	std::vector<int> blockTags;
	std::vector<std::pair<int, int>> refineSeeds;

};

#endif