#include "PlaneSegment.h"

PlaneSegment::PlaneSegment(const PointCloud& points, const int root_block_id,
	const int row_start, const int col_start,
	const int width, const int height,
	const int windowWidth, const int windowHeight,
	const Parameters& params)
{
	root_id = root_block_id;

	bool valid = true;
	for (int i = row_start, icnt = 0; icnt<windowHeight && i<height; ++i, ++icnt) {
		for (int j = col_start, jcnt = 0; jcnt<windowWidth && j<width; ++j, ++jcnt) {
			double x = 0, y = 0, z = 10000;
			if (!points.get(i, j, x, y, z)) {
				valid = false;
				break;
			}
			double xn = 0, yn = 0, zn = 10000;
			if (j + 1<width && (points.get(i, j + 1, xn, yn, zn)
				&& depthDisContinuous(z, zn, params))) {
				valid = false;
				break;
			}
			if (i + 1<height && (points.get(i + 1, j, xn, yn, zn)
				&& depthDisContinuous(z, zn, params))) {
				valid = false;
				break;
			}
			planeData.insert(x, y, z);
		}
		if (!valid)
			break;
	}
	if (valid) {
		merged = false;
		num = planeData.num;
	}
	else {
		num = 0;
		planeData.zero();
		merged = true;
	}

	if (num<4) {
		mse = curvature = std::numeric_limits<double>::quiet_NaN();
	}
	else {
		planeData.GeneratePlane(center, normal, mse, curvature);
	}
}

PlaneSegment::PlaneSegment(const PlaneSegment& pa, const PlaneSegment& pb) : planeData(pa.planeData, pb.planeData)
{
	merged = false;
	root_id = pa.num >= pb.num ? pa.root_id : pb.root_id;
	num = planeData.num;
	planeData.GeneratePlane(center, normal, mse, curvature);
}

double PlaneSegment::simCosine(const PlaneSegment& p) const {
	return std::abs(normal[0] * p.normal[0] +
		normal[1] * p.normal[1] +
		normal[2] * p.normal[2]);
}

double PlaneSegment::distance(const double pt[3]) const {
	return normal[0] * (pt[0] - center[0]) +
		normal[1] * (pt[1] - center[1]) +
		normal[2] * (pt[2] - center[2]);
}

void PlaneSegment::connect(PlaneSegment::Ptr p) {
	if (p) {
		neighbors.insert(p);
		p->neighbors.insert(this);
	}
}

void PlaneSegment::Split() {
	for (auto itr = neighbors.begin(); itr != neighbors.end(); ++itr) {
		PlaneSegment::Ptr neighbor = (*itr);
		neighbor->neighbors.erase(this);
	}
	neighbors.clear();
}

void PlaneSegment::Merge(PlaneSegment& pa, PlaneSegment& pb, DisjointSet& disjoint_set) {
	disjoint_set.Merge(pa.root_id, pb.root_id);

	neighbors.insert(pa.neighbors.begin(), pa.neighbors.end());
	neighbors.insert(pb.neighbors.begin(), pb.neighbors.end());
	neighbors.erase(&pa);
	neighbors.erase(&pb);

	pa.Split();
	pb.Split();

	for (auto itr = neighbors.begin(); itr != neighbors.end(); ++itr) {
		(*itr)->neighbors.insert(this);
	}

	pa.merged = pb.merged = true;
}
