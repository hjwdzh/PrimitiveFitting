#ifndef PLANE_UTIL_H_
#define PLANE_UTIL_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <opencv2/opencv.hpp>

#define DEGREE2RAD(d) ((d)*3.14159265358979323846/180.0)

inline static bool EIG_SVD(double K[3][3], double s[3], double V[3][3])
{
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(
		Eigen::Map<Eigen::Matrix3d>(K[0], 3, 3));
	Eigen::Map<Eigen::Vector3d>(s, 3, 1) = es.eigenvalues();
	Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(V[0], 3, 3) = es.eigenvectors();
	return true;
}

inline static bool EIG_SVD2(double K[2][2], double s[2], double V[2][2])
{
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(
		Eigen::Map<Eigen::Matrix2d>(K[0], 2, 2));
	Eigen::Map<Eigen::Vector2d>(s, 2, 1) = es.eigenvalues();
	Eigen::Map<Eigen::Matrix<double, 2, 2, Eigen::RowMajor>>(V[0], 2, 2) = es.eigenvectors();
	return true;
}

class Parameters {
public:
	enum Phase { INIT = 0, MERGE = 1, REFINE = 2 };

	double d_sigma;
	double epsilon_init;
	double epsilon_merge;

	double z_min, z_max;
	double angle_min, angle_max;
	double sim_merge;
	double sim_refine;

	double alphaTwo;
	double depthChange;

	Parameters() : d_sigma(1.6e-6),
		epsilon_init(5), epsilon_merge(8),
		z_min(500), z_max(4000),
		angle_min(DEGREE2RAD(15.0)), angle_max(DEGREE2RAD(90.0)),
		sim_merge(std::cos(DEGREE2RAD(60.0))),
		sim_refine(std::cos(DEGREE2RAD(15.0))),
		alphaTwo(0.04), depthChange(0.02)
	{}

	inline double ThMSE(const Phase phase, const double z = 0) const {
		switch (phase) {
		case INIT:
			return std::pow(d_sigma*z*z + epsilon_init, 2);
		case MERGE:
		case REFINE:
		default:
			return std::pow(d_sigma*z*z + epsilon_merge, 2);
		}
	}

	inline double ThAng(const Phase phase, const double z = 0) const {
		switch (phase) {
		case INIT:
		{
			double clipped_z = z;
			clipped_z = std::max(clipped_z, z_min);
			clipped_z = std::min(clipped_z, z_max);
			const double factor = (angle_max - angle_min) / (z_max - z_min);
			return std::cos(factor*clipped_z + angle_min - factor*z_min);
		}
		case MERGE:
		{
			return sim_merge;
		}
		case REFINE:
		default:
		{
			return sim_refine;
		}
		}
	}

	inline double ThDz(const double z) const {
		return alphaTwo * fabs(z) + depthChange;
	}
};

inline static bool depthDisContinuous(const double d0, const double d1, const Parameters& params)
{
	return fabs(d0 - d1) > params.ThDz(d0);
}

#include <vector>

class DisjointSet
{
private:
	std::vector<int> parent;
	std::vector<int> rank;

public:
	DisjointSet(const int n)
	{
		parent.resize(n);
		rank.resize(n);

		for (int i = 0; i < n; ++i) {
			parent[i] = i;
			rank[i] = 1;
		}
	}

	~DisjointSet() {}

	inline void remove(const int x) {
		if (parent[x] != x) {
			--rank[GetParent(x)];
			parent[x] = x;
			rank[x] = 1;
		}
	}

	inline int Rank(const int x) {
		return rank[GetParent(x)];
	}

	inline int Merge(const int x, const int y)
	{
		const int xr = GetParent(x);
		const int yr = GetParent(y);

		if (xr == yr)
			return xr;

		const int xrSize = rank[xr];
		const int yrSize = rank[yr];

		if (xrSize < yrSize) {
			parent[xr] = yr;
			rank[yr] += rank[xr];
			return yr;
		}
		else {
			parent[yr] = xr;
			rank[xr] += rank[yr];
			return xr;
		}
	}

	inline int GetParent(const int x)
	{
		if (parent[x] != x)
			parent[x] = GetParent(parent[x]);

		return parent[x];
	}
};

class PointCloud {
public:
	cv::Mat depth;

	PointCloud(const cv::Mat& d, double fx, double fy, double cx, double cy)
	{
		m_fx = fx;
		m_fy = fy;
		m_cx = cx;
		m_cy = cy;
		depth = d;
	}
	int width() const { return depth.cols; }
	int height() const { return depth.rows; }
	bool get(const int row, const int col, double& x, double& y, double& z) const {
		float d = depth.at<float>(row, col);
		if (d == 0)
			return 0;
		z = d;
		x = (col - m_cx) / m_fx * z;
		y = (row - m_cy) / m_fy * z;
		return true;
	}
	double m_fx, m_fy, m_cx, m_cy;
};
																							  
#endif
