#include "PlaneData.h"
#include "PlaneUtil.h"

PlaneData::PlaneData()
{
	zero();
}

PlaneData::PlaneData(const PlaneData& a, const PlaneData& b) :
	sum_x(a.sum_x + b.sum_x), sum_y(a.sum_y + b.sum_y), sum_z(a.sum_z + b.sum_z),
	sum_xx(a.sum_xx + b.sum_xx), sum_yy(a.sum_yy + b.sum_yy), sum_zz(a.sum_zz + b.sum_zz),
	sum_xy(a.sum_xy + b.sum_xy), sum_yz(a.sum_yz + b.sum_yz), sum_xz(a.sum_xz + b.sum_xz), num(a.num + b.num)
{}

void PlaneData::zero() {
	sum_x = sum_y = sum_z = sum_xx = sum_yy = sum_zz = sum_xy = sum_yz = sum_xz = 0;
	num = 0;
}

void PlaneData::insert(const double x, const double y, const double z) {
	sum_x += x; sum_y += y; sum_z += z;
	sum_xx += x*x; sum_yy += y*y; sum_zz += z*z;
	sum_xy += x*y; sum_yz += y*z; sum_xz += x*z;
	++num;
}

void PlaneData::GeneratePlane(double center[3], double normal[3],
	double& mse, double& curvature) const
{
	const double sc = ((double)1.0) / this->num;

	center[0] = sum_x*sc;
	center[1] = sum_y*sc;
	center[2] = sum_z*sc;
	double K[3][3] = {
		{ sum_xx - sum_x*sum_x*sc, sum_xy - sum_x*sum_y*sc, sum_xz - sum_x*sum_z*sc },
		{ 0, sum_yy - sum_y*sum_y*sc, sum_yz - sum_y*sum_z*sc },
		{ 0, 0, sum_zz - sum_z*sum_z*sc }
	};
	K[1][0] = K[0][1]; K[2][0] = K[0][2]; K[2][1] = K[1][2];
	double sv[3] = { 0, 0, 0 };
	double V[3][3] = { 0 };
	EIG_SVD(K, sv, V);

	if (V[0][0] * center[0] + V[1][0] * center[1] + V[2][0] * center[2] <= 0) {
		normal[0] = V[0][0];
		normal[1] = V[1][0];
		normal[2] = V[2][0];
	}
	else {
		normal[0] = -V[0][0];
		normal[1] = -V[1][0];
		normal[2] = -V[2][0];
	}
	mse = sv[0] * sc;
	curvature = sv[0] / (sv[0] + sv[1] + sv[2]);
}

