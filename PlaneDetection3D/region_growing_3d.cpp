/*
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2030. All rights reserved.
 * Description:  Implementation of different 3D Plane Detection Algorithm.
 * Author: Created by huangjingwei 589411
 * Create date: 2021-08-13
 */
// STL includes.
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Random.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <CGAL/Timer.h>

#include <Eigen/Core>
#include <Eigen/Dense>

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using Input_range = CGAL::Point_set_3<Point_3>;
using Point_map = typename Input_range::Point_map;
using Normal_map = typename Input_range::Vector_map;
using Neighbor_query =
    CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, Input_range,
                                                       Point_map>;
// using Neighbor_query =
// CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, Input_range,
// Point_map>;
using Region_type =
    CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<
        Kernel, Input_range, Point_map, Normal_map>;
using Region_growing =
    CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query,
                                          Region_type>;
using Indices = std::vector<std::size_t>;
using Output_range = CGAL::Point_set_3<Point_3>;
using Points_3 = std::vector<Point_3>;

typedef std::pair<Point_3, Vector_3> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixD;
typedef Eigen::Matrix<FT, 3, 1> Vector3;

struct Pointcloud {
  MatrixD P, N, C;
};

void ComputePointNormals(
    Pointcloud& pc,                         // input points + output normals
    unsigned int nb_neighbors_pca_normals)  // number of neighbors
{
  PointList points(pc.P.rows());
  for (int i = 0; i < points.size(); ++i) {
    auto p = pc.P.row(i);
    points[i].first = Point_3(p[0], p[1], p[2]);
  }
  CGAL::Timer task_timer;
  task_timer.start();

  // Estimates normals direction.
  // Note: pca_estimate_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  CGAL::pca_estimate_normals<Concurrency_tag>(
      points, nb_neighbors_pca_normals,
      CGAL::parameters::point_map(
          CGAL::First_of_pair_property_map<PointVectorPair>())
          .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  pc.N.conservativeResize(points.size(), 3);
  for (int i = 0; i < points.size(); ++i) {
    pc.P.row(i) =
        Vector3(points[i].first.x(), points[i].first.y(), points[i].first.z());
    pc.N.row(i) = Vector3(points[i].second.x(), points[i].second.y(),
                          points[i].second.z());
    if (pc.N.row(i).dot(Vector3(1, 1, 1)) < 0) pc.N.row(i) = -pc.N.row(i);
  }
}

// Define an insert iterator. 
int PlaneDetectRegion(Pointcloud& pc,
                      std::vector<std::pair<Vector3, FT> >& plane_parameters,
                      std::vector<int>& new_instances,
                      double dist_thres, double angle_thres, int min_points,
                      int num_neigbhors) {
  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);
  for (int i = 0; i < pc.P.rows(); ++i) {
    input_range.insert(Kernel::Point_3(pc.P(i, 0), pc.P(i, 1), pc.P(i, 2)));
  }
  auto it = input_range.begin();
  for (int i = 0; i < pc.N.rows(); ++i) {
    input_range.normal(*(it++)) =
        Kernel::Vector_3(pc.N(i, 0), pc.N(i, 1), pc.N(i, 2));
  }
  // Default parameter values for the data file point_set_3.xyz.
  const std::size_t k = num_neigbhors;
  const FT max_distance_to_plane = dist_thres;
  const FT max_accepted_angle = angle_thres;
  const std::size_t min_region_size = min_points;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(input_range, k, input_range.point_map());
  Region_type region_type(input_range, max_distance_to_plane,
                          max_accepted_angle, min_region_size,
                          input_range.point_map(), input_range.normal_map());
  // Create an instance of the region growing class.
  Region_growing region_growing(input_range, neighbor_query, region_type);
  // Run the algorithm.
  Output_range output_range;
  std::size_t number_of_regions = 0;

  std::vector<std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  plane_parameters.resize(regions.size());
  new_instances.resize(pc.P.rows(), -1);

  for (int i = 0; i < regions.size(); ++i) {
    Vector3 c(0, 0, 0);
    for (auto& idx : regions[i]) {
      new_instances[idx] = i;
      c += pc.P.row(idx);
    }
    c /= (double)regions[i].size();

    MatrixD C = MatrixD::Zero(3, 3);
    for (auto& idx : regions[i]) {
      Vector3 diff = pc.P.row(idx);
      diff -= c;
      C += diff * diff.transpose();
    }
    Eigen::JacobiSVD<MatrixD> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector3 N = svd.matrixU().col(2);
    plane_parameters[i] = std::make_pair(N, -N.dot(c));
  }
  return regions.size();
}

int main(int argc, char** argv) {

  const int kNeighbors = 16;
  const double kDistThres = 1e-1;
  const double kAngleThres = 20;
  const int kMinPoints = 200;

  // simulate a box
  std::vector<Vector3> points;
  for (int i = 0; i < 100; ++i) {
    for (int j = 0; j < 100; ++j) {
      points.push_back(Vector3(i * 0.01, j * 0.01, 0));
      points.push_back(Vector3(i * 0.01, j * 0.01, 1));
      points.push_back(Vector3(i * 0.01, 0, j * 0.01));
      points.push_back(Vector3(i * 0.01, 1, j * 0.01));
      points.push_back(Vector3(0, i * 0.01, j * 0.01));
      points.push_back(Vector3(1, i * 0.01, j * 0.01));
    }
  }
  for (auto& p : points) {
    double dx = rand() % 256 / 256.0 - 0.5;
    double dy = rand() % 256 / 256.0 - 0.5;
    double dz = rand() % 256 / 256.0 - 0.5;
    p += Vector3(dx * 0.01, dy * 0.01, dz * 0.01);
  }

  Pointcloud pc;
  pc.P.resize(points.size(), 3);
  memcpy(pc.P.data(), points.data(), sizeof(Vector3) * points.size());

  ComputePointNormals(pc, kNeighbors);

  std::vector<std::pair<Vector3, double> > plane_parameters;
  std::vector<int> plane_instances;
  int num_inst = PlaneDetectRegion(pc, plane_parameters, plane_instances,
    kDistThres, kAngleThres, kMinPoints, kNeighbors);

  printf("Num inst: %d\n", num_inst);
  printf("Params:\n");
  for (auto& p : plane_parameters) {
    printf("%f %f %f %f\n", p.first[0], p.first[1], p.first[2], p.second);
  }
  std::vector<Vector3> colors(num_inst);
  for (auto& c : colors) {
    c = Vector3(rand() % 256 / 256.0,
      rand() % 256 / 256.0, rand() % 256 / 256.0);
  }

  std::ofstream os("result.obj");
  for (int i = 0; i < points.size(); ++i) {
    auto p = points[i];
    auto c = Vector3(0, 0, 0);
    if (plane_instances[i] >= 0)
      c = colors[plane_instances[i]];
    os << "v " << p[0] << " " << p[1] << " " << p[2] << " "
      << c[0] << " " << c[1] << " " << c[2] << "\n";
  }
  os.close();
  return 0;
}