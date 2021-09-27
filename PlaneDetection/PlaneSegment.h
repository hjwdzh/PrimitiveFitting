#ifndef PLANE_SEG_H_
#define PLANE_SEG_H_

#include <set>
#include <vector>
#include <limits>

#include <memory>
#include "PlaneUtil.h"
#include "PlaneData.h"

using std::shared_ptr;

class PlaneSegment {
public:
	typedef PlaneSegment* Ptr;
	typedef std::shared_ptr<PlaneSegment> shared_ptr;

	PlaneData planeData;

	std::set<PlaneSegment*> neighbors;

	PlaneSegment(const PointCloud& points, const int root_block_id,
		const int row_start, const int col_start,
		const int width, const int height,
		const int windowWidth, const int windowHeight,
		const Parameters& params);

	PlaneSegment(const PlaneSegment& pa, const PlaneSegment& pb);

	double simCosine(const PlaneSegment& p) const;

	double distance(const double pt[3]) const;

	void connect(PlaneSegment::Ptr p);

	void Split();

	void Merge(PlaneSegment& pa, PlaneSegment& pb, DisjointSet& disjoint_set);

	int root_id;
	double mse;
	double center[3];
	double normal[3];
	double curvature;
	int num;
	bool merged;
};

class SizeHelper {
public:
	bool operator()(const PlaneSegment::shared_ptr& a,
		const PlaneSegment::shared_ptr& b) const {
		return b->num < a->num;
	}
};

class MSEHelper {
public:
	bool operator()(const PlaneSegment::shared_ptr& a,
		const PlaneSegment::shared_ptr& b) const {
		return b->mse < a->mse;
	}
};

#endif
