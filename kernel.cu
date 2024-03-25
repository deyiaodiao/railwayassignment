//The following code implemented by Oleg Konings in association with Morgan Hough and Gazzaley lab
//A simple implementation of the Floyd-Warshall all-pairs-shortest path algorithm with path reconstruction. This is indended to be used on directed graphs with no negative cycles
//The Adjacency Matrix is in Row-major format, and is implemented both in CUDA on a Nvidia GTX 680 2GB GPU, and in serial CPU code using an Intel i7-3770 3.9 ghz.
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <cuda.h>
#include <ctime>
#include <cassert>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define pb push_back 
#define all(c) (c).begin(),(c).end()
#include <Windows.h>
#include <MMSystem.h>
#pragma comment(lib, "winmm.lib")
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>//to detect host memory leaks
using namespace std;


#define _DTH cudaMemcpyDeviceToHost
#define _HTD cudaMemcpyHostToDevice

//these can be altered on user depending on data set and type of operation(random test, read from file etc)
#define BLOCK_SIZE 256
#define RANGE 997
#define RANDOM_GSIZE 4
#define FILE_GSIZE 8298//the number of edges in Wiki-Vote.txt if the file test is run
#define INF (1<<30)
#define DO_TEST_RANDOM 1
#define DO_TEST_FROM_FILE 0

//typedef for vector used in path reconstruction
typedef pair<pair<int, int>, int> Piii;

//forward function declarations
bool InitMMTimer(UINT wTimerRes);
void DestroyMMTimer(UINT wTimerRes, bool init);
void _CPU_Floyd(int* G, int* Gpath, int N);
void _showPath(int start, int end, const vector<Piii>& path, const int* D, const int N);
bool _getPath(int curEdge, int nxtEdge, vector<Piii>& path, const int* D, const int* Dpath, const int N);
void _get_full_paths(const int* D, const int* Dpath, const int N);

//CUDA GPU kernel/functions forward declaration
__global__ void _Wake_GPU(int reps);
__global__ void _GPU_Floyd_kernel(int k, int* G, int* P, int N);
void _GPU_Floyd(int* H_G, int* H_Gpath, const int N);

//other optional utility functions
int _read_from_file(int* G, const int N);
void _generateRandomGraph(int* G, int N, int range, int density);
void _generate_result_file(bool success, unsigned int cpu_time, unsigned int gpu_time, int N);

//------------------------------------------------------------------------------------------------------------------------------------------------------------


// railwayassignment.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include "load_example.h"
#include <queue>
#include <optional>
#include <iomanip>
#include <cmath>
#include "graph.h"
#include "resultnbpr.h"
#include <unordered_set>  
#include <set>  
#include <chrono>
#include <sstream>
#include <random>
#include <omp.h>
#include <algorithm> // 包含 std::find

////网络节点总数
//int N = 5000;
////网络线路链总数
//int M = 800;

////网络节点总数
//int N = 20;
////网络线路链总数
//int M = 4;

//网络节点总数
int N = 1053;
//网络线路链总数
int M = 500;

//更新概率
const int update_rate = 5;

//写入矩阵至csv
const int writefile = 0;

//读取矩阵自csv
const int readfile =0;
//收敛精度
const double ConvergenceAccuracy = 0.00000001;
//最大迭代次数
const int MaxIte = 100;

const double step_size = 1;

//关键节点车站集合
set<int> KeyNodes;
//网络线路链集合
vector<vector<int>> chains(M);
//网络缩减线路链集合
vector<vector<int>> reducedchains(M);
// 全局变量，存储节点
std::vector<Node> Nodes;
//网络映射关系
vector<int> graphMap;

Graph generatedGraph(N);

unordered_map<string, pair<vector<int>, double>> pathsAndWeights;

//------------------------------------------------------------------------------------------------------------------------------------------------------------
//int main() {
//	char ch;
//	srand(time(NULL));
//
//	if (DO_TEST_RANDOM) {//will use the #define(s) to init a random adjacency Matrix of RANDOM_GSIZE size
//		const int NumBytes = RANDOM_GSIZE * RANDOM_GSIZE * sizeof(int);
//		//host allocations to create Adjancency matrix and result matrices with path matrices
//		int* OrigGraph = (int*)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
//		int* H_G = (int*)malloc(NumBytes);
//		int* H_Gpath = (int*)malloc(NumBytes);
//		int* D_G = (int*)malloc(NumBytes);
//		int* D_Gpath = (int*)malloc(NumBytes);
//
//		_generateRandomGraph(OrigGraph, RANDOM_GSIZE, RANGE, 25);//init graph with values
//
//		cout << "Successfully created random highly connected graph in adjacency Matrix form with " << RANDOM_GSIZE * RANDOM_GSIZE << " elements.\n";
//		cout << "Also created 2 pairs of distinct result Matrices to store the respective results of the CPU results and the GPU results.\n";
//		for (int i = 0; i < RANDOM_GSIZE * RANDOM_GSIZE; i++) {//copy for use in computation
//			H_G[i] = D_G[i] = OrigGraph[i];//copy for use in computation
//			H_Gpath[i] = D_Gpath[i] = -1;//set to all negative ones for use in path construction
//		}
//		unsigned int cpu_time = 0, gpu_time = 0;
//		cout << "\nFloyd-Warshall on CPU underway:\n";
//		UINT wTimerRes = 0;
//		bool init = InitMMTimer(wTimerRes);
//		DWORD startTime = timeGetTime();
//
//		//_CPU_Floyd(H_G, H_Gpath, RANDOM_GSIZE);//find shortest paths (with path construction) on serial CPU (Intel i7 3770 3.9 ghz)
//
//		DWORD endTime = timeGetTime();
//		cpu_time = unsigned int(endTime - startTime);
//		printf("CPU Timing: %dms\n", cpu_time);
//		DestroyMMTimer(wTimerRes, init);
//		//wake up GPU from idle
//		cout << "\nFloyd-Warshall on GPU underway:\n";
//		_Wake_GPU << <1, BLOCK_SIZE >> > (32);
//
//		//call host function which will copy all info to device and run CUDA kernels
//		wTimerRes = 0;
//		init = InitMMTimer(wTimerRes);
//		startTime = timeGetTime();
//
//		_GPU_Floyd(D_G, D_Gpath, RANDOM_GSIZE);
//
//		endTime = timeGetTime();
//		gpu_time = unsigned int(endTime - startTime);
//		printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
//		DestroyMMTimer(wTimerRes, init);
//
//		//compare the device generated result against the host generated result
//		cout << "Verifying results of final adjacency Matrix and Path Matrix.\n";
//
//		int same_adj_Matrix = memcmp(H_G, D_G, NumBytes);
//		if (same_adj_Matrix == 0) {
//			cout << "Adjacency Matrices Equal!\n";
//		}
//		else
//			cout << "Adjacency Matrices Not Equal!\n";
//
//		int same_path_Matrix = memcmp(H_Gpath, D_Gpath, NumBytes);
//		if (same_path_Matrix == 0) {
//			cout << "Path reconstruction Matrices Equal!\n";
//		}
//		else
//			cout << "Path reconstruction Matrices Not Equal!\n";
//
//		_get_full_paths(D_G, D_Gpath, RANDOM_GSIZE);//find out exact step-by-step shortest paths between vertices(if such a path exists)
//
//		_generate_result_file(bool(same_adj_Matrix == 0 && same_path_Matrix == 0), cpu_time, gpu_time, RANDOM_GSIZE);
//
//		free(OrigGraph);
//		free(H_G);
//		free(H_Gpath);
//		free(D_G);
//		free(D_Gpath);
//	}
//
//	_CrtDumpMemoryLeaks();
//	cin >> ch;
//	return 0;
//}

bool InitMMTimer(UINT wTimerRes) {
	TIMECAPS tc;
	if (timeGetDevCaps(&tc, sizeof(TIMECAPS)) != TIMERR_NOERROR) { return false; }
	wTimerRes = min(max(tc.wPeriodMin, 1), tc.wPeriodMax);
	timeBeginPeriod(wTimerRes);
	return true;
}

void DestroyMMTimer(UINT wTimerRes, bool init) {
	if (init)
		timeEndPeriod(wTimerRes);
}

void _CPU_Floyd(int* G, int* Gpath, int N) {//standard N^3 algo
	for (int k = 0; k < N; ++k)for (int i = 0; i < N; ++i)for (int j = 0; j < N; ++j) {
		int curloc = i * N + j, loca = i * N + k, locb = k * N + j;
		if (G[curloc] > (G[loca] + G[locb])) {
			G[curloc] = (G[loca] + G[locb]);
			Gpath[curloc] = k;
		}
	}
}

void _showPath(int start, int end, const vector<Piii>& path, const int* D, const int N) {
	cout << "\nHere is the shortest cost path from " << start << " to " << end << ", at a total cost of " << D[start * N + end] << ".\n";
	for (int i = path.size() - 1; i >= 0; --i) {
		cout << "From " << path[i].first.first << " to " << path[i].first.second << " at a cost of " << path[i].second << '\n';
	}
	cout << '\n';
}

bool _getPath(int curEdge, int nxtEdge, vector<Piii>& path, const int* D, const int* Dpath, const int N) {
	int curIdx = curEdge * N + nxtEdge;
	if (D[curIdx] >= INF)return false;
	if (Dpath[curIdx] == -1) {//end of backwards retracement
		path.push_back(make_pair(make_pair(curEdge, nxtEdge), D[curIdx]));
		return true;
	}
	else {//record last edge cost and move backwards
		//path.push_back(make_pair(make_pair(Dpath[curIdx], nxtEdge), D[Dpath[curIdx] * N + nxtEdge]));

		path.push_back(make_pair(make_pair(curEdge,Dpath[curIdx] ), D[curEdge * N + Dpath[curIdx]]));

		//return _getPath(curEdge, Dpath[curIdx], path, D, Dpath, N);
		return _getPath(Dpath[curIdx], nxtEdge, path, D, Dpath, N);
	}
}

void _get_full_paths(const int* D, const int* Dpath, const int N) {
	int start_vertex = -1, end_vertex = -1;
	vector<Piii> path;
	do {
		path.clear();
		cout << "Enter start vertex #:";
		cin >> start_vertex;
		cout << "Enter dest vertex(enter negative number to exit) #:";
		cin >> end_vertex;
		if (start_vertex < 0 || start_vertex >= N || end_vertex < 0 || end_vertex >= N)return;

		if (_getPath(start_vertex, end_vertex, path, D, Dpath, N)) {
			_showPath(start_vertex, end_vertex, path, D, N);

		}
		else {
			cout << "\nThere does not exist valid a path between " << start_vertex << " , and " << end_vertex << '\n';

		}
	} while (1);
}

__global__ void _Wake_GPU(int reps) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= reps)return;
}

__global__ void _GPU_Floyd_kernel(int k, int *G,int *P, int N){//G will be the adjacency matrix, P will be path matrix
	int col=blockIdx.x*blockDim.x + threadIdx.x;
	if(col>=N)return;
	int idx=N*blockIdx.y+col;

	__shared__ int best;
	if(threadIdx.x==0)
		best=G[N*blockIdx.y+k];
	__syncthreads();
	if(best==INF)return;
	int tmp_b=G[k*N+col];
	if(tmp_b==INF)return;
	int cur=best+tmp_b;
	if(cur<G[idx]){
		G[idx]=cur;
		P[idx]=k;
	}
}

void _GPU_Floyd(int* H_G, int* H_Gpath, const int N) {
	cout << "_GPU_Floyd!\n";
	//allocate device memory and copy graph data from host
	int* dG, * dP;
	int numBytes = N * N * sizeof(int);
	cudaError_t err = cudaMalloc((int**)&dG, numBytes);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
	err = cudaMalloc((int**)&dP, numBytes);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
	//copy from host to device graph info
	err = cudaMemcpy(dG, H_G, numBytes, _HTD);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
	err = cudaMemcpy(dP, H_Gpath, numBytes, _HTD);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }

	dim3 dimGrid((N + BLOCK_SIZE - 1) / BLOCK_SIZE, N);

	for (int k = 0; k < N; k++) {//main loop

		_GPU_Floyd_kernel << <dimGrid, BLOCK_SIZE >> > (k, dG, dP, N);
		err = cudaThreadSynchronize();
		if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
	}
	//copy back memory
	err = cudaMemcpy(H_G, dG, numBytes, _DTH);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
	err = cudaMemcpy(H_Gpath, dP, numBytes, _DTH);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }

	//free device memory
	err = cudaFree(dG);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
	err = cudaFree(dP);
	if (err != cudaSuccess) { printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); }
}

void _generateRandomGraph(int* G, int N, int range, int density) {//density will be between 0 and 100, indication the % of number of directed edges in graph
	//range will be the range of edge weighting of directed edges
	//int Prange = (100 / density);

	int Prange = 100;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {//set G[i][i]=0
				G[i * N + j] = 0;
				continue;
			}
			//int pr = rand() % Prange;
			//G[i * N + j] = pr == 0 ? ((rand() % range) + 1) : INF;//set edge random edge weight to random value, or to INF

			int pr = rand() % Prange+1;
			G[i * N + j] = pr;//set edge random edge weight to random value
		}
	}
}

int _read_from_file(int* G, const int N) {//reads in edge list from file
	int num_edges = 0;

	ifstream readfile;//enable stream for reading file
	readfile.open("Wiki-Vote.txt");
	assert(readfile.good());//make sure it finds the file & file is
	string line;
	int v0, v1;
	while (getline(readfile, line)) {
		istringstream linestream(line);
		linestream >> v0 >> v1;
		G[v0 * N + v1] = 1;
		num_edges++;
	}
	readfile.close();
	return num_edges;
}

void _generate_result_file(bool success, unsigned int cpu_time, unsigned int gpu_time, int N) {

	if (!success) {
		cout << "Error in calculation!\n";
		return;
	}
	else {
		ofstream myfile;
		myfile.open("Floyd-Warshall_result.txt");
		myfile << "Success! The GPU Floyd-Warshall result and the CPU Floyd-Warshall results are identical(both final adjacency matrix and path matrix).\n\n";
		myfile << "N= " << N << " , and the total number of elements(for Adjacency Matrix and Path Matrix) was " << N * N << " .\n";
		myfile << "Matrices are int full dense format(row major) with a minimum of " << (N * N) / 4 << " valid directed edges.\n\n";
		myfile << "The CPU timing for all was " << float(cpu_time) / 1000.0f << " seconds, and the GPU timing(including all device memory operations(allocations,copies etc) ) for all was " << float(gpu_time) / 1000.0f << " seconds.\n";
		myfile << "The GPU result was " << float(cpu_time) / float(gpu_time) << " faster than the CPU version.\n";
		myfile.close();
	}
}


//----------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool isElementInSet(const std::set<T>& mySet, const T& element) {
	// 使用 find 函数查找元素
	auto it = mySet.find(element);
	// 判断元素是否存在于 set 中
	return (it != mySet.end());
}

double generateRandomdouble(double min, double max) {
	random_device rd;  // 获取一个随机设备种子
	mt19937 gen(rd()); // 使用Mersenne Twister算法生成随机数引擎
	uniform_real_distribution<double> dis(min, max); // 使用均匀分布生成随机浮点数

	return dis(gen);
}


 //打印邻接矩阵
void printAdjacencyMatrix(const vector<vector<int>>& graph) {
	for (const auto& row : graph) {
		for (int node : row) {
			cout << node << " ";
		}
		cout << endl;
	}
}



void printPath(vector<vector<int>>& next, int u, int v) {
	if (next[u][v] == -1) {
		cout << "No path exists between " << u << " and " << v << endl;
		return;
	}
	if (u == 0 && v == 2479) {
		cout << endl;
	}
	cout << "Shortest path from " << u << " to " << v << " is: " << u << " ";
	while (u != v) {
		u = next[u][v];
		//cout << "-> " << u << " "<<"next " << next[u][v] << " ";

		cout << "-> " << u << " ";
	}
	cout << endl;
}

std::vector<std::vector<std::vector<int>>> getShortestPaths(const std::vector<std::vector<int>>& next) {
	int numNodes = next.size();
	std::vector<std::vector<std::vector<int>>> paths(numNodes, std::vector<std::vector<int>>(numNodes));

	for (int i = 0; i < numNodes; ++i) {
		for (int j = 0; j < numNodes; ++j) {
			if (i != j && next[i][j] != -1) {
				int source = i;
				int destination = j;
				std::vector<int> path;
				while (source != destination) {
					path.push_back(source);
					source = next[source][destination];
				}
				path.push_back(destination);
				paths[i][j] = path;
			}
		}
	}

	return paths;
}

void printPaths(const std::vector<std::vector<std::vector<int>>>& paths) {
	int numNodes = paths.size();
	for (int i = 0; i < numNodes; ++i) {
		for (int j = 0; j < numNodes; ++j) {
			if (paths[i][j].empty()) {
				std::cout << "No path exists from " << i << " to " << j << std::endl;
			}
			else {
				std::cout << "Shortest path from " << i << " to " << j << ": ";
				for (int node : paths[i][j]) {
					std::cout << node << " -> ";
				}
				std::cout << j << std::endl;
			}
		}
	}
}

bool hasNegativeCycle(const std::vector<std::vector<int>>& next) {
	int numNodes = next.size();

	for (int i = 0; i < numNodes; ++i) {
		if (next[i][i] != -1) {
			cout << i << "存在负权回路 " << endl;
			return true; // 存在负权回路
		}
	}
	return false; // 不存在负权回路
}

vector<vector<int>> floydWarshall(vector<vector<double>>& adjMatrix, int numNodes) {
	//omp_set_num_threads(12); // 设置使用线程	
	// // 获取当前时间点
	auto start_time = std::chrono::high_resolution_clock::now();


	cout << "Floyd-Warshall begin " << endl;
	// 初始化距离矩阵和next矩阵
	vector<vector<double>> dist(numNodes, vector<double>(numNodes));
	vector<vector<int>> next(numNodes, vector<int>(numNodes, -1));

#pragma omp parallel for
	for (int i = 0; i < numNodes; ++i) {
		for (int j = 0; j < numNodes; ++j) {
			dist[i][j] = adjMatrix[i][j];
			if (adjMatrix[i][j] != INF && i != j) {
				next[i][j] = j;
			}
		}
	}

	// 开始计算最短路径
//#pragma omp parallel for
	for (int k = 0; k < numNodes; ++k) {
		//std::cout << "Update: " << i << std::endl;
		//if (k % 100 == 0) {
		//	std::cout << "FLoyd: " << k << std::endl;
		//}

#pragma omp parallel for
		for (int i = 0; i < numNodes; ++i) {
			//if (i % 100 == 0) {
			//	std::cout << "FLoyd: " << k << "    "<<i << std::endl;
			//}
			for (int j = 0; j < numNodes; ++j) {
				if (dist[i][k] == INF || dist[k][j] == INF) continue;
#pragma omp critical
				if (dist[i][k] + dist[k][j] < dist[i][j]) {
					dist[i][j] = dist[i][k] + dist[k][j];
					next[i][j] = next[i][k]; // 更新next矩阵
				}
				//else if (dist[i][k] + dist[k][j] == dist[i][j] && next[i][k] != -1) {
				//	next[i][j] = next[i][k]; // 检查是否有更短的路径可用
				//}
			}
		}
	}
	cout << "Floyd-Warshall over" << endl;
	// 输出最短路径
	//cout << "最短路径距离矩阵：" << endl;
	//for (int i = 0; i < numNodes; ++i) {
	//	for (int j = 0; j < numNodes; ++j) {
	//		if (dist[i][j] == INF)
	//			cout << "INF\t";
	//		else
	//			cout << dist[i][j] << "\t";
	//	}
	//	cout << endl;
	//}

	//std::vector<std::vector<std::vector<int>>> paths = getShortestPaths(next);

	//if (hasNegativeCycle(next)) {
	//	cout << "存在负权回路 " << endl;
	//}

	 //输出最短路径节点序列
	//cout << endl << "最短路径节点序列：" << endl;
	//for (int i = 0; i < numNodes; ++i) {

	//	cout << "Shortest path from " << i << " to " << endl;
	//	for (int j = 0; j < numNodes; ++j) {
	//		if (i != j&& next[i][j] != -1) {
	//			//cout << "Shortest path from " << i << " to " << j << ": ";
	//			printPath(next, i, j);
	//		}
	//	}
	//}

	// 获取当前时间点

	auto end_time = std::chrono::high_resolution_clock::now();

	// 计算时
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// 将时间转换为秒，并保留三位小数
	double durationInSeconds = duration.count() / 1000.0;

	// 输出执行时间
	std::cout << "FW程序执行时间: " << std::fixed << std::setprecision(3) << durationInSeconds << "秒" << std::endl;






	return next;
}

// 写入邻接矩阵到 CSV 文件
void writeAdjMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& adjMatrix) {
	std::ofstream file(filename + ".csv");
	if (file.is_open()) {
		for (const auto& row : adjMatrix) {
			for (size_t i = 0; i < row.size(); ++i) {
				file << row[i];
				if (i != row.size() - 1) {
					file << ",";
				}
			}
			file << std::endl;
		}
		file.close();
		std::cout << "Adjacency matrix has been written to " << filename << ".csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file " << filename << ".csv" << " for writing." << std::endl;
	}
}


void writeNextMatrixToCSV(const std::string& filename, const std::vector<std::vector<int>>& _Next) {
	std::ofstream file(filename + ".csv");
	if (file.is_open()) {
		for (const auto& row : _Next) {
			for (size_t i = 0; i < row.size(); ++i) {
				file << row[i];
				if (i != row.size() - 1) {
					file << ",";
				}
			}
			file << std::endl;
		}
		file.close();
		std::cout << "Adjacency matrix has been written to " << filename << ".csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file " << filename << ".csv" << " for writing." << std::endl;
	}
}


// 从 CSV 文件中读取邻接矩阵
std::vector<std::vector<double>> readAdjMatrixFromCSV(const std::string& filename) {
	std::vector<std::vector<double>> adjMatrix;
	std::ifstream file(filename + ".csv");
	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {
			std::vector<double> row;
			std::stringstream ss(line);
			std::string cell;
			while (std::getline(ss, cell, ',')) {
				row.push_back(std::stod(cell));
			}
			adjMatrix.push_back(row);
		}
		file.close();
		std::cout << "Adjacency matrix has been read from " << filename << ".csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file " << filename << ".csv" << " for reading." << std::endl;
	}
	return adjMatrix;
}

std::vector<std::vector<Path>> getAllPaths(const std::vector<std::vector<int>>& next) {
	int numNodes = next.size();
	std::vector<std::vector<Path>> allPaths(numNodes, std::vector<Path>(numNodes));

	for (int i = 0; i < numNodes; ++i) {
		for (int j = 0; j < numNodes; ++j) {
			if (i != j && next[i][j] != -1) {
				int source = i;
				int destination = j;
				std::vector<int> path;
				while (source != destination) {
					path.push_back(source);
					source = next[source][destination];
				}
				path.push_back(destination);

				// 存储路径到结构体 Path 中
				Path p;
				p.path = path;
				allPaths[i][j] = p;
			}
		}
	}

	return allPaths;
}

Optional<Path> getShortestPath(const std::vector<std::vector<int>>& next, int src, int dest) {
	std::vector<int> path;
	int source = src;
	int destination = dest;

	while (source != destination) {
		//std::cout << "--> " << destination << std::endl;
		path.push_back(source);
		source = next[source][destination];
		if (source == -1) {
			// 如果没有路径连接源节点和目标节点，则返回空
			return Optional<Path>();
		}
	}
	path.push_back(destination);

	// 存储路径到结构体 Path 中
	Path p;
	p.path = path;

	return Optional<Path>(p);
}


Optional<Path> getShortestPathGPU(const std::vector<std::vector<int>>& next,
	std::vector<std::vector<double>>& _adjMatrix, 
	int src, int dest) {
	std::vector<int> path;

	std::vector<int> reverse_path;
	int source = src;
	int destination = dest;

	while (source != destination) {
		if (_adjMatrix[source][next[source][destination]] < INF) {
			path.push_back(source);
			source = next[source][destination];
		}
		//std::cout << "--> " << destination << std::endl;
		if (_adjMatrix[next[source][destination]][destination] < INF) {
			reverse_path.push_back(destination);
			destination = next[source][destination];
		}
		if (source == -1) {
			// 如果没有路径连接源节点和目标节点，则返回空
			return Optional<Path>();
		}
		if (next[source][destination] == destination) {
			path.push_back(source);
			reverse_path.push_back(destination);
			break;
			// 如果没有路径连接源节点和目标节点，则返回空
			//return Optional<Path>();
		}
		if (_adjMatrix[next[source][destination]][destination] >= INF
			&& _adjMatrix[source][next[source][destination]] >=INF) {
			Optional<Path> oppath1= getShortestPathGPU(next, _adjMatrix, source, next[source][destination]);
			Optional<Path> oppath2 = getShortestPathGPU(next, _adjMatrix, next[source][destination], destination);
			if (oppath1.has_value() && oppath2.has_value())
				// 修改 ResultODPaths 的 paths 和 flow
			{
				std::vector<int> path1 = oppath1.get_value().path;
				std::vector<int> path2 = oppath2.get_value().path;
				path.insert(path.end(), path1.begin(), path1.end());
				path.insert(path.end(), path2.begin() + 1, path2.end());
			}
			else {
				return Optional<Path>();
			}
			break;
			// 如果没有路径连接源节点和目标节点，则返回空
			//return Optional<Path>();
		}
	}

	std::reverse(reverse_path.begin(), reverse_path.end());

	// 将第二个vector添加到第一个vector末尾
	path.insert(path.end(), reverse_path.begin(), reverse_path.end());

	for (size_t i = 0; i < path.size()-1; ++i) {
		if (_adjMatrix[path[i]][path[i + 1]] >= INF) {
			cout << path[i] << "\t" << path[i + 1] << "\tINF link\t" << src <<"\t"<<dest << endl;
		}
	}


	// 存储路径到结构体 Path 中
	Path p;
	p.path = path;

	return Optional<Path>(p);
}



void computeAllOrNothing(std::vector<std::vector<ODPaths>>& ResultODPaths,
	const std::vector<std::vector<int>>& ReducedNext, 
	std::vector<std::vector<double>>& _adjMatrix,
	int N) {
	// 遍历所有源点和目标点
	for (int i = 0; i < N; ++i) {
		if (i % 100 == 0) {
			std::cout << "AON: " << i << std::endl;
		}
		//std::cout << "AON: " << i << std::endl;
		for (int j = 0; j < N; ++j) {
			// 计算最短路径
			if (i != j && ReducedNext[i][j] != -1) {
				//Optional<Path> shortestPath = getShortestPath(ReducedNext, i, j);
				//if (i == 2 && j == 59) {
				//	int dafn = 0; 
				//}
				Optional<Path> shortestPath = getShortestPathGPU(ReducedNext, _adjMatrix, i, j);
				if (shortestPath.has_value())
				// 修改 ResultODPaths 的 paths 和 flow
				{
					Path path = shortestPath.get_value();
					ResultODPaths[i][j].paths.push_back(path);
					ResultODPaths[i][j].flow.push_back(ResultODPaths[i][j].demand);
				}
			}
		}
	}
}


//需保证矩阵中最大值为INF
vector<vector<int>> _gpu_fw_apsp(std::vector<std::vector<double>>& _adjMatrix,int _N) {



	// 获取当前时间点
	auto start_time = std::chrono::high_resolution_clock::now();


	vector<vector<int>>_Next;

	const int NumBytes = N * N * sizeof(int);
	int* D_G = (int*)malloc(NumBytes);
	int* D_Gpath = (int*)malloc(NumBytes);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if(_adjMatrix[i][j]!= INF){ 
				D_G[i * N + j] = _adjMatrix[i][j] * 1000; 
				D_Gpath[i * N + j] = j;
			}
			else { 
				D_G[i * N + j] = INF; 
				D_Gpath[i * N + j] = -1;	
			}
			//D_Gpath[i * N + j] = -1;
			// 
			//D_G[i * N + j] = _adjMatrix[i][j];
			
		}
	}
	_Wake_GPU << <1, BLOCK_SIZE >> > (32);
	_GPU_Floyd(D_G, D_Gpath, N);
	_Next.clear(); // 清空 _Next

	// 为 _Next 分配空间
	_Next.resize(N, std::vector<int>(N));

	// 将 D_Gpath 中的值复制到 _Next 中
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			_Next[i][j] = D_Gpath[i * N + j];
		}
	}
	
	//_get_full_paths(D_G, D_Gpath, N);//find out exact step-by-step shortest paths between vertices(if such a path exists)


	free(D_G);
	free(D_Gpath);

	// 获取当前时间点

	auto end_time = std::chrono::high_resolution_clock::now();

	// 计算时
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// 将时间转换为秒，并保留三位小数
	double durationInSeconds = duration.count() / 1000.0;

	// 输出执行时间
	std::cout << "GPUFW时间: " << std::fixed << std::setprecision(3) << durationInSeconds << "秒" << std::endl;


	return _Next;
}







int main()
{
	//vector<vector<double>> matrix01 = readAdjMatrixFromCSV(to_string(N) + "GPUMLBPR01");
	//vector<vector<double>> matrix02 = readAdjMatrixFromCSV(to_string(N) + "GPUMLBPR02");
	//vector<vector<double>> matrix03 = readAdjMatrixFromCSV(to_string(N) + "GPUMLBPR03");
	//vector<vector<double>> matrix04 = readAdjMatrixFromCSV(to_string(N) + "GPUMLBPR04");
	//vector<vector<double>> matrix05 = readAdjMatrixFromCSV(to_string(N) + "GPUMLBPR05");

	//vector<vector<double>> matrix06 = readAdjMatrixFromCSV(to_string(N) + "GPUMLNEXT01");
	//vector<vector<double>> matrix07 = readAdjMatrixFromCSV(to_string(N) + "GPUMLNEXT02");
	//vector<vector<double>> matrix08 = readAdjMatrixFromCSV(to_string(N) + "GPUMLNEXT03");
	//vector<vector<double>> matrix09 = readAdjMatrixFromCSV(to_string(N) + "GPUMLNEXT04");
	//vector<vector<double>> matrix10 = readAdjMatrixFromCSV(to_string(N) + "GPUMLNEXT05");
	//for (int i = 0; i < N; ++i) {
	//	for (int j = 0; j < N; ++j) {
	//		//cout << i<<"\t" << j <<endl;
	//		if (std::fabs(matrix02[i][j] - matrix01[i][j]) > 1e-19) { std::cout << i << "2\t" << j << "\t" << matrix02[i][j] << "\t" << matrix01[i][j] << std::endl; }
	//		if (std::fabs(matrix03[i][j] - matrix01[i][j]) > 1e-19) { std::cout << i << "3\t" << j << "\t" << matrix03[i][j] << "\t" << matrix01[i][j] << std::endl; }
	//		if (std::fabs(matrix04[i][j] - matrix01[i][j]) > 1e-19) { std::cout << i << "4\t" << j << "\t" << matrix04[i][j] << "\t" << matrix01[i][j] << std::endl; }
	//		if (std::fabs(matrix05[i][j] - matrix01[i][j]) > 1e-19) { std::cout << i << "5\t" << j << "\t" << matrix05[i][j] << "\t" << matrix01[i][j] << std::endl; }
	//		if (std::fabs(matrix07[i][j] - matrix06[i][j]) > 1e-19) { std::cout << i << "7\t" << j << "\t" << matrix07[i][j] << "\t" << matrix06[i][j] << std::endl; }
	//		if (std::fabs(matrix08[i][j] - matrix06[i][j]) > 1e-19) { std::cout << i << "8\t" << j << "\t" << matrix08[i][j] << "\t" << matrix06[i][j] << std::endl; }
	//		if (std::fabs(matrix09[i][j] - matrix06[i][j]) > 1e-19) { std::cout << i << "9\t" << j << "\t" << matrix09[i][j] << "\t" << matrix06[i][j] << std::endl; }
	//		if (std::fabs(matrix10[i][j] - matrix06[i][j]) > 1e-19) { std::cout << i << "10\t" << j << "\t" << matrix10[i][j] << "\t" << matrix06[i][j] << std::endl; }
	//	}
	//}

	//const int NumBytes = N * N * sizeof(int);
	////host allocations to create Adjancency matrix and result matrices with path matrices
	//int* OrigGraph = (int*)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
	//int* H_G = (int*)malloc(NumBytes);
	//int* H_Gpath = (int*)malloc(NumBytes);
	//int* D_G = (int*)malloc(NumBytes);
	//int* D_Gpath = (int*)malloc(NumBytes);
	//_generateRandomGraph(OrigGraph, N, RANGE, 25);//init graph with values
	//for (int i = 0; i < N * N; i++) {//copy for use in computation
	//	H_G[i] = D_G[i] = OrigGraph[i];//copy for use in computation
	//	H_Gpath[i] = D_Gpath[i] = -1;//set to all negative ones for use in path construction
	//}

	//std::string filename = "Philadelphia_trips - 副本.csv"; // Replace with your CSV file path
	//Philadelphia_trips
	//int num_demands = 1525; // Example value for num_demands
	//int num_net_nodes = 13390;
	//std::vector<std::vector<double>> data = readPHITripDataFromCSV(filename, num_demands);

	std::string filename = "Winnipeg_trips - 副本.csv"; // Replace with your CSV file path
	//Winnipeg_trips
	int num_demands = 147; // Example value for num_demands
	//int num_net_nodes = 1053;

	std::vector<std::vector<double>> Winnipeg_trips_data = readWinTripDataFromCSV(filename, num_demands);

	int zeros = 0;
	int non_zeros = 0;
	for (int i = 1; i < 2; ++i) {
		for (int j = 0; j < num_demands + 1; ++j) {
			// cout << i << " " << j << " " << data[i][j] << endl;
			if (Winnipeg_trips_data[i][j] == 0) {
				zeros++;
			}
			else {
				non_zeros++;
			}
		}
	}
	cout << zeros << " " << non_zeros << endl;


	//std::string file_path = "Philadelphia_net - 副本.csv"; // 请将文件路径替换为实际文件路径
	std::string file_path = "Winnipeg_net - 副本.csv"; // 请将文件路径替换为实际文件路径

	std::vector<Link> links = read_data(file_path);


	cout << "the number of nodes (N): " << N << endl;

	cout << "Initializing begin!" << endl;

	Graph ReducedGraph(N);
	//ReducedGraph.generateRandomWeights();
	// 
	// 遍历 links 向量
	for (const auto& link : links) {
		ReducedGraph.adjMatrix[link.init_node][link.term_node] = link.free_flow_time;
		//processLink(link);  // 处理当前 Link 对象
	}



	if (writefile == 1) {

		writeAdjMatrixToCSV(to_string(N) + to_string(M) + "规模", ReducedGraph.adjMatrix);
	}
	if (readfile == 1) {
		ReducedGraph.adjMatrix = readAdjMatrixFromCSV(to_string(N) + to_string(M) + "规模");
	}

	////cpu实现fw最短路
	//vector<vector<int>> ReducedNext_check = floydWarshall(ReducedGraph.adjMatrix, N);
	vector<vector<int>> ReducedNext=_gpu_fw_apsp(ReducedGraph.adjMatrix, N);

	//vector<vector<int>> ReducedNext1 = _gpu_fw_apsp(ReducedGraph.adjMatrix, N);
	
	

	//for (int i = 0; i < N; ++i) {
	//	for (int j = 0; j < N; ++j) {
	//		//cout << i<<"\t" << j <<endl;
	//		if (std::fabs(ReducedNext[i][j] - ReducedNext1[i][j]) > 1e-19) { std::cout << i << "\t" << j << "\t" << ReducedNext[i][j] << "\t" << ReducedNext[i][j] << std::endl; }
	//		//cout << i << "\t" << j << "\t" << ReducedNext[i][j] << "\t" << ReducedNext[i][j] << std::endl;
	//	}
	//}

	//for (int i = 0; i < N; ++i) {
	//	for (int j = 0; j < N; ++j) {
	//		cout << "ID\t " << i << "\t" << j << '\n';
	//		cout << "cost\t " << ReducedGraph.adjMatrix[i][j] << "\t" << D_G[i * N + j] << '\n';
	//		//cout <<   "Next\t " << ReducedNext_check[i][j] << "\tD_Gpath\t" << D_Gpath[i * N + j] << endl;

	//		cout << "Next\t " << ReducedNext[i][j] << '\t' << ReducedNext_check[i][j] << "\tD_Gpath\t" << D_Gpath[i * N + j] << endl;
	//	}
	//}






	// 获取当前时间点
	auto start_time = std::chrono::high_resolution_clock::now();



	cout << "Initializing get ODPaths!" << endl;

	//GraphResult ReducedResult(N, N, ReducedGraph.adjMatrix);
	//std::vector<std::vector<ODPaths>> ResultODPaths(N, std::vector<ODPaths>(N, ODPaths(3)));

	GraphResult ReducedResult(N, N, ReducedGraph.adjMatrix);
	std::vector<std::vector<ODPaths>> ResultODPaths(N, std::vector<ODPaths>(N, ODPaths(0)));
	for (int i = 0; i < num_demands; ++i) {
		for (int j = 0; j < num_demands; ++j) {
			ResultODPaths[i][j].demand = Winnipeg_trips_data[i][j];
		}
	}






	cout << "Initializing all or nothing!" << endl;
	//1 初始化
	//1.1 初始化全有全无
	computeAllOrNothing(ResultODPaths, ReducedNext, ReducedGraph.adjMatrix, N);

	cout << "Initializing flow!" << endl;
	//1.2 根据全有全无结果更新flow、bpr、diff_bpr
	for (int i = 0; i < N; ++i) {
		//std::cout << "Update: " << i << std::endl;
		if (i % 1000 == 0) {
			std::cout << "Update: " << i << std::endl;
		}
		for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
			if (i != j) {
				//update_flow_bpr_derivative(ResultODPaths[i][j],
				//	ReducedResult.flow,
				//	ReducedResult.diff_bpr,
				//	ReducedResult.bpr,
				//	ReducedResult.alpha,
				//	ReducedResult.beta,
				//	ReducedResult.capacity,
				//	ReducedGraph.adjMatrix);
				//更新flow
				for (int k = 0; k < ResultODPaths[i][j].paths.size(); ++k) {
					const Path& path = ResultODPaths[i][j].paths[k];
					double singel_flow = ResultODPaths[i][j].flow[k];
					for (int t = 0; t < path.path.size() - 1; ++t) {
						int from = path.path[t];
						int to = path.path[t + 1];
						ReducedResult.flow[from][to] += singel_flow;
					}
				}
				//更新bpr、diff_bpr
				for (int k = 0; k < ResultODPaths[i][j].paths.size(); ++k) {
					const Path& path = ResultODPaths[i][j].paths[k];
					double singel_flow = ResultODPaths[i][j].flow[k];
					for (int t = 0; t < path.path.size() - 1; ++t) {
						int from = path.path[t];
						int to = path.path[t + 1];
						ReducedResult.bpr[from][to] = bpr_function(ReducedGraph.adjMatrix[from][to],
							ReducedResult.alpha[from][to],
							ReducedResult.beta[from][to],
							ReducedResult.flow[from][to],
							ReducedResult.capacity[from][to]);
						ReducedResult.diff_bpr[from][to] = bpr_function_derivative(ReducedGraph.adjMatrix[from][to],
							ReducedResult.alpha[from][to],
							ReducedResult.beta[from][to],
							ReducedResult.flow[from][to],
							ReducedResult.capacity[from][to]);
					}
				}
			}
		}
	}
//
//
	//2 主循环
	int current_iteration = 0; // 当前迭代次数
	while (current_iteration < MaxIte) {
		std::cout << "Iteration: " << current_iteration << "-------------------------------------------------------------------------" << std::endl;
		//if (current_iteration % 100 == 0) {
		//	std::cout << "Iteration: " << current_iteration << std::endl;
		//}
		//vector<vector<int>> loopReducedNext = floydWarshall(ReducedResult.bpr, N);
		vector<vector<int>> loopReducedNext = _gpu_fw_apsp(ReducedResult.bpr, N);
		//vector<vector<int>> loopReducedNext = _gpu_fw_apsp(ReducedGraph.adjMatrix, N);






		std::cout << "Iteration: " << current_iteration << " each OD " << std::endl;
//#pragma omp parallel for
				//2.1 对每个OD求解子问题
		for (int i = 0; i < N; ++i) {
			//if (i % 100 == 0) {
			//	std::cout << "Main loop i: " << i << std::endl;
			//}
			//std::cout << "Main loop i: " << i << std::endl;

			int remainder = (i+ current_iteration) % update_rate;
			if (remainder != 0) { continue; }

			for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
				//if (j % 100 == 0) {
				//	std::cout << "Main loop j: " << j << std::endl;
				//}
				//std::cout << "Main loop j: " << j << std::endl;
				if (i != j) {
					if (ResultODPaths[i][j].demand == 0) { continue; }
					////根据概率决定是否更新
					//std::random_device rd;
					//// 使用随机设备生成引擎
					//std::mt19937 gen(rd());
					//// 定义 double 分布，范围为0到1
					//std::uniform_real_distribution<double> dis(0.0, 1.0);

					//// 生成随机 double 值并打印
					//double random_double = dis(gen);
					//if (random_double > update_rate) { continue; }





					//2.1.1 列生成
					//判断路径是否存在于 paths 中，不存在则添加到末尾
					//if (i == 2 && j == 59) {
					//	std::cout << "Iteration: " << current_iteration << " Adjustment " << i << "   " << j << "   " << std::endl;
					//}

					Optional<Path> loopSP = getShortestPath(loopReducedNext, i, j);
					if (loopSP.has_value()) {
						//std::cout << "Optional contains a value." << std::endl;
						add_path_if_not_exists(ResultODPaths[i][j], loopSP.get_value());
					}
					if (ResultODPaths[i][j].paths.size() == 1) { continue; }
					//2.1.2 基于当前列集合重新分配流量
					//2.1.2.1 选定标准路径
					ResultODPaths[i][j].path_bpr.clear(); // 清空原有的内容
					//计算所有路径bpr
					for (const Path& path : ResultODPaths[i][j].paths) {
						double path_bpr = calculate_total_bpr(path.path, ReducedResult.bpr); // 获得当前路径对应的 BPR 值
						ResultODPaths[i][j].path_bpr.push_back(path_bpr);

					}
					std::pair<double, int> min_value_and_index = get_min_value_and_index(ResultODPaths[i][j].path_bpr);
					double min_bpr = min_value_and_index.first;
					int min_index = min_value_and_index.second;

					//2.1.2.2 针对每条路径，计算与标准路径之间的优化方案，并更新各路径流量
					//计算一阶导
					ResultODPaths[i][j].path_diff_bpr.clear(); // 清空原有的内容
					for (const Path& path : ResultODPaths[i][j].paths) {
						double path_diff_bpr = calculate_total_diff_bpr(path.path, ResultODPaths[i][j].paths[min_index].path, ReducedResult.diff_bpr); // 获得当前路径对应的 BPR 值
						ResultODPaths[i][j].path_diff_bpr.push_back(path_diff_bpr);

					}

					//计算流量调整值 步长*当前路径流量*（当前路径总阻抗函数值-基本路径总阻抗函数值）/当前路径总阻抗一阶导
					for (size_t path_id = 0; path_id < ResultODPaths[i][j].paths.size(); ++path_id) {
						if (ResultODPaths[i][j].path_diff_bpr[path_id] == 0) { continue; }
						if (path_id == min_index) { continue; }
						double Adjustment = step_size *
							ResultODPaths[i][j].flow[path_id] *
							(ResultODPaths[i][j].path_bpr[path_id] - ResultODPaths[i][j].path_bpr[min_index]) /
							ResultODPaths[i][j].path_diff_bpr[path_id];
						//std::cout << "Iteration: " << current_iteration << " Adjustment " <<i<<"   "<<j<<"   "<< Adjustment << std::endl;
						if (Adjustment > ResultODPaths[i][j].flow[path_id]) {
							Adjustment = ResultODPaths[i][j].flow[path_id];
						}
						//std::cout << "Iteration: " << current_iteration << " Adjustment " << i << "   " << j << "   " << Adjustment << std::endl;
						ResultODPaths[i][j].flow[path_id] = ResultODPaths[i][j].flow[path_id] - Adjustment;
						ResultODPaths[i][j].flow[min_index] = ResultODPaths[i][j].flow[min_index] + Adjustment;
					}
					//2.1.3 判定单个OD是否满足收敛性指标（可选）
				}
			}
		}

		std::cout << "Iteration: " << current_iteration << " update flow " << std::endl;
		//根据调整结果，更新flow、bpr、diff_bpr    
#pragma omp parallel for
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
				if (i != j) {
					ReducedResult.flow[i][j] = 0.0f;
					ReducedResult.bpr[i][j] = 0.0f;
					ReducedResult.diff_bpr[i][j] = 0.0f;
				}
			}
		}

		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
				if (i != j) {
					//更新flow
					for (int k = 0; k < ResultODPaths[i][j].paths.size(); ++k) {
						const Path& path = ResultODPaths[i][j].paths[k];
						double singel_flow = ResultODPaths[i][j].flow[k];
						for (int t = 0; t < path.path.size() - 1; ++t) {
							int from = path.path[t];
							int to = path.path[t + 1];
							ReducedResult.flow[from][to] += singel_flow;
						}
					}
					//更新bpr、diff_bpr
					for (int k = 0; k < ResultODPaths[i][j].paths.size(); ++k) {
						const Path& path = ResultODPaths[i][j].paths[k];
						double singel_flow = ResultODPaths[i][j].flow[k];
						for (int t = 0; t < path.path.size() - 1; ++t) {
							int from = path.path[t];
							int to = path.path[t + 1];
							ReducedResult.bpr[from][to] = bpr_function(ReducedGraph.adjMatrix[from][to],
								ReducedResult.alpha[from][to],
								ReducedResult.beta[from][to],
								ReducedResult.flow[from][to],
								ReducedResult.capacity[from][to]);
							ReducedResult.diff_bpr[from][to] = bpr_function_derivative(ReducedGraph.adjMatrix[from][to],
								ReducedResult.alpha[from][to],
								ReducedResult.beta[from][to],
								ReducedResult.flow[from][to],
								ReducedResult.capacity[from][to]);
						}
					}
				}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
				if (i != j) {
					if (ReducedResult.bpr[i][j] == 0.0f) { ReducedResult.bpr[i][j] = INF; }
				}
			}
		}

		double total_convergence = calculate_total_convergence(ResultODPaths, ReducedResult.bpr);
		std::cout << std::fixed; // 设置固定的小数点表示
		std::cout << std::setprecision(10); // 设置小数点后20位
		std::cout << "total_convergence: " << total_convergence << std::endl;

		//for (int i = 0; i < N; ++i) {
		//	for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
		//		if (ResultODPaths[i][j].demand==0) { continue; }
		//		if (i != j) {
		//			std::cout << "i : " << i<<" j "<<j<<" "<< ResultODPaths[i][j].demand << std::endl;
		//			for (int k = 0; k < ResultODPaths[i][j].paths.size();++k) {
		//				std::cout << "flow: " << ResultODPaths[i][j].flow[k] <<"  "<< ResultODPaths[i][j].path_bpr[k] << std::endl;
		//			}

		//		}
		//	}
		//}

		if (total_convergence > 1 - ConvergenceAccuracy) { break; }
		// 增加当前迭代次数
		current_iteration++;








		// 可以在这里添加其他终止条件的判断逻辑
	}


//	//


	// 获取当前时间点

	auto end_time = std::chrono::high_resolution_clock::now();

	// 计算时
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// 将时间转换为秒，并保留三位小数
	double durationInSeconds = duration.count() / 1000.0;

	// 输出执行时间
	std::cout << "程序执行时间: " << std::fixed << std::setprecision(3) << durationInSeconds << "秒" << std::endl;



	//free(OrigGraph);
	//free(H_G);
	//free(H_Gpath);
	//free(D_G);
	//free(D_Gpath);






	//printAdjacencyMatrix(generatedGraph);
	system("pause");
	return 0;
}

