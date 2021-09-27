#ifndef PLANE_DATA_H_
#define PLANE_DATA_H_

class PlaneData {
public:
	PlaneData();
	PlaneData(const PlaneData& a, const PlaneData& b);
	
	void zero();
	void insert(const double x, const double y, const double z);

	void GeneratePlane(double center[3], double normal[3], double& mse, double& curvature) const;

	double sum_x, sum_y, sum_z, sum_xx, sum_yy, sum_zz, sum_xy, sum_yz, sum_xz;
	int num;

};
#endif