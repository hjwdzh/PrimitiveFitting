#include "PlaneFitter.h"


PlaneFitter::PlaneFitter() : 
	points(0), width(0), height(0),
	max_group_iter(100000), min_block_size(3000),
	windowWidth(10), windowHeight(10)
{
}

PlaneFitter::~PlaneFitter() {}

void PlaneFitter::clear() {
	points = 0;
	planes.clear();
	disjoint_set.reset();
	root_block.clear();
	blockTags.clear();
	refineSeeds.clear();
}

int PlaneFitter::Neighbors(const int i, const int j, const int H, const int W, int neighbors[4])
{
	const int id = i*W + j;
	int cnt = 0;
	if (j>0) neighbors[cnt++] = (id - 1);
	if (j<W - 1) neighbors[cnt++] = (id + 1);
	if (i>0) neighbors[cnt++] = (id - W);
	if (i<H - 1) neighbors[cnt++] = (id + W);
	return cnt;
}


int PlaneFitter::BlockID(const int cx, const int cy) const {
//	assert(pixX >= 0 && pixY >= 0 && pixX<width && pixY<height);
	const int block_w = width / windowWidth;
	const int block_h = height / windowHeight;
	const int by = cy / windowHeight;
	const int bx = cx / windowWidth;
	return (by<block_h && bx<block_w) ? (by*block_w + bx) : -1;
}

void PlaneFitter::InitGraph(QueueMSE& queueMSE) {
	const int block_h = height / windowHeight;
	const int block_w = width / windowWidth;

	std::vector<PlaneSegment::Ptr> graph(block_h * block_w, 0);

	for (int i = 0; i < block_h; ++i) {
		for (int j = 0; j < block_w; ++j) {
			PlaneSegment::shared_ptr p(new PlaneSegment(
				*points, i * block_w + j,
				i * windowHeight, j * windowWidth,
				width, height,
				windowWidth, windowHeight,
				params));
			if (p->mse < params.ThMSE(Parameters::INIT, p->center[2])
				&& !p->merged)
			{
				graph[i * block_w + j] = p.get();
				queueMSE.push(p);
			}
			else {
				graph[i * block_w + j] = 0;
			}
		}
	}

	for (int i = 0; i < block_h; ++i) {
		for (int j = 1; j < block_w; j += 2) {
			const int idx = i * block_w + j;
			if (graph[idx - 1] == 0) {
				--j;
				continue;
			}
			if (graph[idx] == 0)
				continue;
			if (j < block_w - 1 && graph[idx + 1] == 0) {
				++j;
				continue;
			}

			const double th_ang = params.ThAng(Parameters::INIT, graph[idx]->center[2]);
			if ((j < block_w - 1 && graph[idx - 1]->simCosine(*graph[idx + 1]) >= th_ang)
				|| (j == block_w - 1 && graph[idx]->simCosine(*graph[idx - 1]) >= th_ang))
			{
				graph[idx]->connect(graph[idx - 1]);
				if (j < block_w - 1) graph[idx]->connect(graph[idx + 1]);
			}
			else {
				--j;
			}
		}
	}

	for (int j = 0; j < block_w; ++j) {
		for (int i = 1; i < block_h; i += 2) {
			const int idx = i * block_w + j;
			if (graph[idx - block_w] == 0) {
				--i;
				continue;
			}
			if (graph[idx] == 0)
				continue;
			if (i < block_h - 1 && graph[idx + block_w] == 0) {
				++i;
				continue;
			}

			const double th_ang = params.ThAng(Parameters::INIT, graph[idx]->center[2]);
			if ((i < block_h - 1 && graph[idx - block_w]->simCosine(*graph[idx + block_w]) >= th_ang)
				|| (i == block_h - 1 && graph[idx]->simCosine(*graph[idx - block_w]) >= th_ang))
			{
				graph[idx]->connect(graph[idx - block_w]);
				if (i<block_h - 1) graph[idx]->connect(graph[idx + block_w]);
			}
			else {
				--i;
			}
		}
	}
}

void PlaneFitter::MergeGraph(QueueMSE& queueMSE) {
	int step = 0;
	while (!queueMSE.empty() && step <= max_group_iter) {
		if (step > max_group_iter)
			break;
		PlaneSegment::shared_ptr p = queueMSE.top();
		queueMSE.pop();
		if (p->merged) {
			assert(p->neighbors.size() <= 0);
			continue;
		}
		PlaneSegment::shared_ptr potential_merge;
		PlaneSegment::Ptr potential_neighbor(0);
		for (auto itr = p->neighbors.begin(); itr != p->neighbors.end(); itr++) {
			PlaneSegment::Ptr neighbor = (*itr);
			if (p->simCosine(*neighbor) < params.ThAng(Parameters::MERGE, p->center[2]))
				continue;
			PlaneSegment::shared_ptr merge(new PlaneSegment(*p, *neighbor));
			if (potential_merge == 0 || potential_merge->mse > merge->mse ||
				(potential_merge->mse == merge->mse && potential_merge->num < merge->num))
			{
				potential_merge = merge;
				potential_neighbor = neighbor;
			}
		}
		if (potential_merge != 0 && potential_merge->mse<params.ThMSE(
			Parameters::MERGE, potential_merge->center[2]))
		{
			queueMSE.push(potential_merge);
			potential_merge->Merge(*p, *potential_neighbor, *disjoint_set);
		}
		else {
			if (p->num >= min_block_size) {
				planes.push_back(p);
			}
			p->Split();
		}
		++step;
	}
	while (!queueMSE.empty()) {
		const PlaneSegment::shared_ptr p = queueMSE.top();
		queueMSE.pop();
		if (p->num >= min_block_size) {
			planes.push_back(p);
		}
		p->Split();
	}
	static SizeHelper sizecmp;
	std::sort(planes.begin(),
		planes.end(),
		sizecmp);
}

void PlaneFitter::AnalyzeBlockTags(std::vector<bool>& validPlanes) {
	root_block.clear();
	for (int block_tag = 0; block_tag<(int)planes.size(); ++block_tag) {
		root_block.insert(std::pair<int, int>(planes[block_tag]->root_id, block_tag));
	}

	const int block_h = height / windowHeight;
	const int block_w = width / windowWidth;
	const int block_size = windowHeight*windowWidth;

	tags.create(height, width, CV_32SC1);
	tags.setTo(-1);
	blockTags.resize(block_h*block_w);

	validPlanes.resize(planes.size(), false);
	int block_id = 0;
	for (int i = 0; i<block_h; ++i) {
		for (int j = 0; j < block_w; ++j) {
			const int setid = disjoint_set->GetParent(block_id);
			const int setSize = disjoint_set->Rank(setid)*block_size;
			if (setSize >= min_block_size) {
				int neighbors[4] = { -1 };
				const int nNbs = Neighbors(i, j, block_h, block_w, neighbors);
				bool flag = true;
				for (int k = 0; k<nNbs; ++k) {
					if (disjoint_set->GetParent(neighbors[k]) != setid)
					{
						flag = false;
						break;
					}
				}
				const int block_tag = root_block[setid];
				if (flag) {
					blockTags[block_id] = block_tag;
					const int by = block_id / block_w;
					const int bx = block_id - by*block_w;
					tags(cv::Range(by*windowHeight, (by + 1)*windowHeight),
						cv::Range(bx*windowWidth, (bx + 1)*windowWidth)).setTo(block_tag);
					validPlanes[block_tag] = true;
				}
				else {
					blockTags[block_id] = -1;
				}
			}
			else {
				blockTags[block_id] = -1;
			}

			if (blockTags[block_id] < 0) {
				if (i > 0) {
					const int upId = block_id - block_w;
					if (blockTags[upId] >= 0) {
						const int upTag = blockTags[upId];
						const int seedCoord = (i*windowHeight - 1)*width + j*windowWidth;
						for (int k = 1; k<windowWidth; ++k) {
							refineSeeds.push_back(std::pair<int, int>(seedCoord + k, upTag));
						}
					}
				}
				if (j > 0) {
					const int leftId = block_id - 1;
					if (blockTags[leftId] >= 0) {
						const int leftTag = blockTags[leftId];
						const int seedCoord = (i*windowHeight)*width + j*windowWidth - 1;
						for (int k = 0; k<windowHeight - 1; ++k) {
							refineSeeds.push_back(std::pair<int, int>(seedCoord + k*width, leftTag));
						}
					}
				}
			}
			else {
				const int block_tag = blockTags[block_id];
				if (i > 0) {
					const int upId = block_id - block_w;
					if (blockTags[upId] != block_tag) {
						const int seedCoord = (i*windowHeight)*width + j*windowWidth;
						for (int k = 0; k<windowWidth - 1; ++k) {
							refineSeeds.push_back(std::pair<int, int>(seedCoord + k, block_tag));
						}
					}
				}
				if (j > 0) {
					const int leftId = block_id - 1;
					if (blockTags[leftId] != block_tag) {
						const int seedCoord = (i*windowHeight)*width + j*windowWidth;
						for (int k = 1; k<windowHeight; ++k) {
							refineSeeds.push_back(std::pair<int, int>(seedCoord + k*width, block_tag));
						}
					}
				}
			}
			++block_id;
		}
	}
}

void PlaneFitter::FloodFill()
{
	std::vector<float> distMap(height*width,
		std::numeric_limits<float>::max());

	for (int k = 0; k < (int)refineSeeds.size(); ++k) {
		const int seedCoord = refineSeeds[k].first;
		const int cy = seedCoord / width;
		const int cx = seedCoord - cy*width;
		const int block_tag = refineSeeds[k].second;
		const PlaneSegment& segment = *planes[block_tag];

		int neighbors[4] = { -1 };
		const int num_neighbors = Neighbors(cy, cx, height, width, neighbors);
		for (int itr = 0; itr < num_neighbors; ++itr) {
			const int idx = neighbors[itr];
			int& trail = tags.at<int>(idx);
			if (trail <= -6) continue;
			if (trail >= 0 && trail == block_tag) continue;
			const int cy = idx / width;
			const int cx = idx - cy * width;
			const int block_id = BlockID(cx, cy);
			if (block_id >= 0 && blockTags[block_id] >= 0) continue;

			double pt[3] = { 0 };
			float cdist = -1;
			if (points->get(cy, cx, pt[0], pt[1], pt[2]) &&
				std::pow(cdist = (float)std::abs(segment.distance(pt)), 2)<9 * segment.mse + 1e-5)
			{
				if (trail >= 0) {
					PlaneSegment& neighbor = *planes[trail];
					if (segment.simCosine(neighbor) >= params.ThAng(Parameters::REFINE, segment.center[2])) {
						neighbor.connect(planes[block_tag].get());
					}
				}
				float& old_dist = distMap[idx];
				if (cdist<old_dist) {
					trail = block_tag;
					old_dist = cdist;
					refineSeeds.push_back(std::pair<int, int>(idx, block_tag));
				}
				else if (trail<0) {
					trail -= 1;
				}
			}
			else {
				if (trail<0) trail -= 1;
			}
		}
	}
}

void PlaneFitter::Process(const PointCloud* _points)
{
	clear();
	points = _points;
	height = points->height();
	width = points->width();
	disjoint_set.reset(new DisjointSet((height / windowHeight) * (width / windowWidth)));

	QueueMSE queueMSE;
	InitGraph(queueMSE);

	MergeGraph(queueMSE);
	Refine();
}

void PlaneFitter::Refine()
{
	std::vector<bool> validPlanes;
	AnalyzeBlockTags(validPlanes);

	std::vector<int> border;

	FloodFill();

	std::vector<PlaneSegment::shared_ptr> originalPlanes;
	planes.swap(originalPlanes);
	QueueMSE queueMSE;
	for (int i = 0; i<(int)originalPlanes.size(); ++i) {
		if (validPlanes[i])
			queueMSE.push(originalPlanes[i]);
	}
	MergeGraph(queueMSE);

	std::vector<int> compactIndices(originalPlanes.size(), -1);
	int offset = 0;
	for (int i = 0; i<(int)originalPlanes.size(); ++i) {
		const PlaneSegment& op = *originalPlanes[i];
		if (!validPlanes[i]) {
			compactIndices[i] = -1;
			continue;
		}
		int np_root_id = disjoint_set->GetParent(op.root_id);
		if (np_root_id == op.root_id) {
			if (compactIndices[i]<0)
				compactIndices[i] = offset++;
		}
		else {
			const int npid = root_block[np_root_id];
			if (compactIndices[npid] < 0) {
				compactIndices[i] = offset;
				compactIndices[npid] = offset++;
			}
			else
				compactIndices[i] = compactIndices[npid];
		}
	}

	static const cv::Vec3b blackColor(0, 0, 0);
	const int nPixels = width * height;
	for (int i = 0; i<nPixels; ++i) {
		int& tag = tags.at<int>(i);
		if (tag >= 0 && compactIndices[tag] >= 0) {
			tag = compactIndices[tag];
		}
	}
}
