#include "Plane.h"

int main(int argc, char** argv) {
	std::string color_path = "../test/color.jpg";
	std::string depth_path = "../test/depth.pgm";
	std::string output_path = "mask.png";

	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], "-c") == 0 && i < argc - 1)
			color_path = argv[i + 1];
		if (strcmp(argv[i], "-d") == 0 && i < argc - 1)
			depth_path = argv[i + 1];
		if (strcmp(argv[i], "-o") == 0) {
			output_path = argv[i + 1];
		}
		if (strcmp(argv[i], "-h") == 0) {
			printf("./test_plane_depth -c color_path -d depth_path -o output_path\n");
		}
	}

	cv::Mat color, depth, mask;
	color = cv::imread(color_path);
	depth = cv::imread(depth_path, cv::IMREAD_UNCHANGED);

	auto res = ExtractPlane(depth, 575.0, 575.0, 319.5, 239.5);
	cv::Mat label = res.first;
	std::vector<PlaneHelper>& plane_param = res.second;

	int num_planes = plane_param.size();
	std::vector<cv::Vec3b> random_color(num_planes);
	for (auto& r : random_color) {
		r = cv::Vec3b(rand() % 256, rand() % 256, rand() % 256);
	}

	for (int i = 0; i < color.rows; ++i) {
		for (int j = 0; j < color.cols; ++j) {
			auto& c = color.at<cv::Vec3b>(i, j);
			int l = label.at<unsigned char>(i, j);
			// 255 means not a plane
			if (l == 255) {
				c = c / 4;
			} else {
				auto plane_color = random_color[l];
				c = c / 4 + plane_color / 4 * 3;
			}
		}
	}

	cv::imwrite(output_path, color);
	return 0;
}
