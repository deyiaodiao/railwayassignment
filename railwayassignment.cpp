// railwayassignment.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>

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
using namespace std;
////网络节点总数
//int N = 5000;
////网络线路链总数
//int M = 800;

////网络节点总数
//int N = 20;
////网络线路链总数
//int M = 4;

//网络节点总数
int N = 1200;
//网络线路链总数
int M = 500;

//写入矩阵至csv
const int writefile = 1;

//读取矩阵自csv
const int readfile =0;
//收敛精度
const double ConvergenceAccuracy =0.00001;
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


// 生成一个包含 n 个顶点的链
vector<int> generateChain(int n,int start) {
	vector<int> chain;
	for (int i = 0; i < n; ++i) {
		//cout << int(i + start);
		chain.push_back(int(i+start));
	}
	return chain;
}

unordered_set<int> getTotalUniqueElements(const std::vector<std::vector<int>>& reducedChains) {
	std::unordered_set<int> uniqueElements;
	for (const auto& chain : reducedChains) {
		for (int element : chain) {
			uniqueElements.insert(element);
		}
	}
	return uniqueElements;
}

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

//缩减线路链规模
vector<int> reduceVector(const vector<int>& original, const set<int>& keyNodes) {
	vector<int> reducedVector;
	// 将第一个关键节点加入新的vector
	reducedVector.push_back(original[0]);
	int pre_key = 0;
	int current_key = -1;

	// 遍历关键节点，处理关键节点之间的中间元素
	for (int i = 1; i < original.size(); ++i) {
		//存在于set中
		if (isElementInSet(keyNodes, original[i])) {
			current_key = i;
			//相邻
			if (current_key - pre_key == 1) {
				reducedVector.push_back(original[i]);
				pre_key = i;
			}
			//不相邻
			else if (current_key - pre_key >= 1)
			{
				int middle_key= round((current_key + pre_key) / 2);
				reducedVector.push_back(original[middle_key]);
				reducedVector.push_back(original[i]);
				pre_key = i;
			}
			else
			{
				std::cout << "缩短链意外情况 " << current_key << "   "<< pre_key << std::endl;
			}
		}
		else {
			continue;
		}

	}

	return reducedVector;
}

// 生成图
Graph generateGraph(int N, int M) {


	srand(time(0));

	// 计算每个链的长度
	int chainLength = N / M;
	int extraNodes = N % M; // 余下的点

	

	for (int i = 0; i < M; ++i) {
		int length = chainLength;
		//if (extraNodes > 0) {
		//	length++;
		//	extraNodes--;
		//}
		chains[i] = generateChain(length, length * i);// 生成一条链
		
		//前后链衔接，保证全联通
		if (i > 0) {
			chains[i].insert(chains[i].begin(), chains[i-1].back());
		}

	}
	for (int i = N - extraNodes; i < N; ++i) {
		chains[0].push_back(i);
	}

	vector<int> KeyNode;
	//添加线路起终点
	for (int i = 0; i < M; ++i) {
		int start = chains[i][0];
		int end = chains[i][chains[i].size() - 1];
		if (i == 0) {
			KeyNode.push_back(start);
			KeyNode.push_back(end);
		}
		auto itstart = std::find(KeyNode.begin(), KeyNode.end(), start);		

		// 检查是否找到了元素
		if (itstart != KeyNode.end()) {
			continue;
			//std::cout << "Vector contains " << target << std::endl;
		}
		else {
			KeyNode.push_back(start);
		}
		auto itend = std::find(KeyNode.begin(), KeyNode.end(), end);
		if (itend != KeyNode.end()) {
			continue;
			//std::cout << "Vector contains " << target << std::endl;
		}
		else {
			KeyNode.push_back(end);
		}


	}


	int numChains30Percent = round(M * 0.3);
	int numChains50Percent = round(M * 0.5);
	int numChains20Percent = round(M * 0.2);
	vector<int> PercentStation30;
	vector<int> PercentStation20;

	//枢纽站 找到30%的链中的一个非首尾元素，插入到另外50%的链中
	for (int i = 0; i < numChains30Percent; ++i) {
		int randomElementIndex = rand() % (chains[i].size() - 2) + 1;
		int selectedElement = chains[i][randomElementIndex];
		// 使用 std::find 查找元素
		auto it = std::find(PercentStation30.begin(), PercentStation30.end(), selectedElement);

		// 找到了元素
		if (it != PercentStation30.end()) {
			continue;
			//std::cout << "Vector contains " << target << std::endl;
		}

		PercentStation30.push_back(selectedElement);



		for (int j = numChains20Percent; j < numChains20Percent + numChains50Percent; ++j) {

			// 使用 std::find 查找元素
			auto it = std::find(chains[j].begin(), chains[j].end(), selectedElement);

			// 检查是否找到了元素
			if (it != chains[j].end()) {
				continue;
				//std::cout << "Vector contains " << target << std::endl;
			}
			else {
				chains[j].push_back(selectedElement);
				//std::cout << "Vector does not contain " << target << std::endl;
			}
		}
	}

	//一般站 找到20%的链中的一个非首尾元素，插入到另外30%的链中
	for (int i = numChains30Percent; i < numChains30Percent+ numChains20Percent; ++i) {
		int randomElementIndex = rand() % (chains[i].size() - 2) + 1;
		int selectedElement = chains[i][randomElementIndex];

		// 使用 std::find 查找元素
		auto it = std::find(PercentStation20.begin(), PercentStation20.end(), selectedElement);

		// 找到了元素
		if (it != PercentStation20.end()) {
			continue;
			//std::cout << "Vector contains " << target << std::endl;
		}


		PercentStation20.push_back(selectedElement);

		for (int j = numChains50Percent; j < numChains50Percent + numChains30Percent; ++j) {

			// 使用 std::find 查找元素
			auto it = std::find(chains[j].begin(), chains[j].end(), selectedElement);

			// 检查是否找到了元素
			if (it != chains[j].end()) {
				continue;
				//std::cout << "Vector contains " << target << std::endl;
			}
			else {
				chains[j].push_back(selectedElement);
				//std::cout << "Vector does not contain " << target << std::endl;
			}
		}
	}


	//新建nodes
	for (int i = 0; i < N; ++i) {
		Node newNode(i);
		Nodes.push_back(newNode);
	}


	// 遍历每个链
	for (int i = 0; i < M; ++i) {
		// 遍历每个链中的每个元素
		for (int j = 0; j < chains[i].size()-1; ++j) {
			generatedGraph.addEdge(chains[i][j], chains[i][j+1], generateRandomdouble(5, 9));
			if (Nodes[chains[i][j]].chain_number == -1) {
				Nodes[chains[i][j]].chain_number = i;
			}
			Nodes[chains[i][j]].edges_A[chains[i][j+1]] = 
			{ chains[i][j + 1], generatedGraph.getWeight(chains[i][j], chains[i][j + 1])};  // 与节点chains[i][j+1]连接
			Nodes[chains[i][j+1]].edges_A[chains[i][j]] =
			{ chains[i][j], generatedGraph.getWeight(chains[i][j+1], chains[i][j]) };  // 与节点chains[i][j]连接

		}
	}

	// 将三个vector的元素插入到set中  
	KeyNodes.insert(KeyNode.begin(), KeyNode.end());
	KeyNodes.insert(PercentStation30.begin(), PercentStation30.end());
	KeyNodes.insert(PercentStation20.begin(), PercentStation20.end());


	return generatedGraph;
}
// 生成缩减图
Graph generateReducedGraph(const std::vector<vector<int>>& chains, const std::set<int>& keyElements) {
	for (int i = 0; i < M; ++i) {
		reducedchains[i] = reduceVector(chains[i], keyElements);
	
	}
	unordered_set<int> totalUnique = getTotalUniqueElements(reducedchains);

	Graph graph_B(totalUnique.size());

	//完成AB网络映射
	int identifier_B = 0;
	graphMap.resize(totalUnique.size());
	std::fill(graphMap.begin(), graphMap.end(), -1);
	for (int element : totalUnique) {
		graphMap[identifier_B] = Nodes[element].identifier;
		Nodes[element].identifier_B = identifier_B;
		++identifier_B;
		//std::cout << element << " "<< Nodes[element].identifier_B << endl;
	}

#pragma omp parallel for
	for (int i = 0; i < totalUnique.size(); ++i) {
		
		for (int j = i+1; j < totalUnique.size(); ++j) {
			std::cout << i << " 当前构造" <<j<< endl;
			//vector<int> path;
			//double weight = generatedGraph.shortestPath(i, j, path);
			//pathsAndWeights[to_string(graphMap[i]) + "-" + to_string(graphMap[j])] = make_pair(path, weight);
			//graph_B.addEdge(i, j, weight);

			graph_B.addEdge(i, j, generateRandomdouble(5, 9));
		}
	}

	return graph_B;
}
// 打印邻接矩阵
//void printAdjacencyMatrix(const vector<vector<int>>& graph) {
//	for (const auto& row : graph) {
//		for (int node : row) {
//			cout << node << " ";
//		}
//		cout << endl;
//	}
//}

// 输出图的dot文件
//void writeDotFile(const vector<vector<int>>& graph, int N, const string& filename) {
//	ofstream file(filename);
//	if (!file.is_open()) {
//		cerr << "Failed to open file for writing." << endl;
//		return;
//	}
//	file << "graph G {" << endl;
//	for (int i = 0; i < N; ++i) {
//		for (int j = i + 1; j < N; ++j) {
//			if (graph[i][j] == 1) {
//				file << i << " -- " << j << ";" << endl;
//			}
//		}
//	}
//	file << "}" << endl;
//	file.close();
//}

// Dijkstra 算法实现
Path dijkstra_shortest_path(const std::vector<std::vector<double>>& adjacency_matrix, int start, int end) {
	int n = adjacency_matrix.size();
	std::vector<double> dist(n, std::numeric_limits<double>::max());
	std::vector<int> prev(n, -1);
	std::vector<bool> visited(n, false);

	dist[start] = 0.0;

	for (int count = 0; count < n - 1; ++count) {
		int u = -1;
		double min_dist = std::numeric_limits<double>::max();

		// 选择未访问的节点中距离最小的节点
		for (int i = 0; i < n; ++i) {
			if (!visited[i] && dist[i] < min_dist) {
				min_dist = dist[i];
				u = i;
			}
		}

		if (u == -1) break;

		visited[u] = true;

		// 更新与节点 u 相邻的节点的距离
		for (int v = 0; v < n; ++v) {
			if (!visited[v] && adjacency_matrix[u][v] && dist[u] + adjacency_matrix[u][v] < dist[v]) {
				dist[v] = dist[u] + adjacency_matrix[u][v];
				prev[v] = u;
			}
		}
	}

	// 从终点向前回溯路径
	Path shortest_path;
	for (int at = end; at != -1; at = prev[at]) {
		shortest_path.path.push_back(at);
	}
	std::reverse(shortest_path.path.begin(), shortest_path.path.end());

	return shortest_path;
}

Path bidirectional_dijkstra(const vector<vector<double>>& adjacency_matrix, int start, int end) {
	int n = adjacency_matrix.size();
	// 从终点向前回溯路径
	Path shortest_path;
	// 初始化正向和反向的距离数组和路径数组
	vector<double> dist_forward(n, INF);
	vector<double> dist_backward(n, INF);
	dist_forward[start] = 0;
	dist_backward[end] = 0;

	// 创建正向和反向的优先队列，保存节点和距离
	priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> forward_queue;
	priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> backward_queue;
	forward_queue.push({ 0, start });
	backward_queue.push({ 0, end });

	// 正向和反向的前驱节点数组
	vector<int> prev_forward(n, -1);
	vector<int> prev_backward(n, -1);

	while (!forward_queue.empty() && !backward_queue.empty()) {
		// 从正向和反向的优先队列中分别取出距离最小的节点
		int u = forward_queue.top().second;
		forward_queue.pop();
		int v = backward_queue.top().second;
		backward_queue.pop();

		// 如果在正向和反向搜索中遇到相同的节点，则返回最短路径
		if (dist_forward[u] + dist_backward[u] < dist_forward[end] || dist_forward[v] + dist_backward[v] < dist_forward[end]) {
			// 回溯路径
			//vector<int> path;
			int node = u;
			while (node != -1) {
				shortest_path.path.push_back(node);
				node = prev_forward[node];
			}
			node = v;
			stack<int> temp_path;
			while (node != -1) {
				temp_path.push(node);
				node = prev_backward[node];
			}
			while (!temp_path.empty()) {
				shortest_path.path.push_back(temp_path.top());
				temp_path.pop();
			}
			return shortest_path;
		}

		// 更新正向和反向的距离数组
		// 正向更新
		for (int i = 0; i < n; ++i) {
			if (dist_forward[i] > dist_forward[u] + adjacency_matrix[u][i]) {
				dist_forward[i] = dist_forward[u] + adjacency_matrix[u][i];
				prev_forward[i] = u;
				forward_queue.push({ dist_forward[i], i });
			}
		}
		// 反向更新
		for (int i = 0; i < n; ++i) {
			if (dist_backward[i] > dist_backward[v] + adjacency_matrix[v][i]) {
				dist_backward[i] = dist_backward[v] + adjacency_matrix[v][i];
				prev_backward[i] = v;
				backward_queue.push({ dist_backward[i], i });
			}
		}
	}

	// 如果无法从起点到达终点，则返回空路径
	return {};
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

		cout << "-> " << u  << " ";
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
			cout <<i<< "存在负权回路 " << endl;
			return true; // 存在负权回路
		}
	}
	return false; // 不存在负权回路
}





vector<vector<int>> floydWarshall(vector<vector<double>>& adjMatrix, int numNodes) {
	//omp_set_num_threads(12); // 设置使用线程	
	// // 获取当前时间点
	auto start_time = std::chrono::high_resolution_clock::now();


	cout << "Floyd-Warshall begin "<<endl;
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

struct Edge {
	int start_node;
	int end_node;
	double capacity;
	double flow;
	double travel_time;
};


std::vector<Edge> initialize_network() {
	// 构造一个简单的5条边的网络
	std::vector<Edge> network = {
		{0, 1, 50.0, 0.0, 1.0},  // 边 0
		{0, 2, 100.0, 0.0, 1.0},  // 边 1
		{1, 3, 100.0, 0.0, 1.0},  // 边 2
		{2, 3, 100.0, 0.0, 1.0},  // 边 3
		{3, 4, 100.0, 0.0, 1.0}   // 边 4
	};
	return network;
}

void user_equilibrium(std::vector<Edge>& network, double t0, double alpha, double beta, double tolerance, int max_iterations) {
	int iteration = 0;
	double delta_flow = 0.0;
	while (iteration < max_iterations) {
		delta_flow = 0.0;
		for (Edge& edge : network) {
			double current_flow = edge.flow;
			double new_flow = edge.capacity / bpr_function(t0, alpha, beta, edge.flow, edge.capacity);
			edge.flow = new_flow;
			delta_flow += std::abs(new_flow - current_flow);
		}
		if (delta_flow < tolerance) {
			break;
		}
		iteration++;
	}
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

std::optional<Path> getShortestPath(const std::vector<std::vector<int>>& next, int src, int dest) {
	std::vector<int> path;
	int source = src;
	int destination = dest;

	while (source != destination) {
		//std::cout << "--> " << destination << std::endl;
		path.push_back(source);
		source = next[source][destination];
		if (source == -1) {
			// 如果没有路径连接源节点和目标节点，则返回空
			return std::nullopt;
		}
	}
	path.push_back(destination);

	// 存储路径到结构体 Path 中
	Path p;
	p.path = path;

	return p;
}


void computeAllOrNothing(std::vector<std::vector<ODPaths>>& ResultODPaths, const std::vector<std::vector<int>>& ReducedNext, int N) {
	// 遍历所有源点和目标点
	for (int i = 0; i < N; ++i) {
		if (i % 10 == 0) {
			std::cout << "AON: " << i << std::endl;
		}
		std::cout << "AON: " << i << std::endl;
		for (int j = 0; j < N; ++j) {
			// 计算最短路径
			if (i != j && ReducedNext[i][j] != -1) {
				std::optional<Path> shortestPath = getShortestPath(ReducedNext, i, j);
				// 修改 ResultODPaths 的 paths 和 flow
				{
					ResultODPaths[i][j].paths.push_back(shortestPath.value());
					ResultODPaths[i][j].flow.push_back(ResultODPaths[i][j].demand);
				}
			}
		}
	}
}





int main()
{





	//int N, M;
	cout << "the number of nodes (N): "<<N<<endl;
	//cin >> N;
	//cout << "the number of chains (M): "<<M<<endl;
	//cin >> M;

	cout << "Initializing begin!"  << endl;

	//N = 500;
	//M = 100;

	
	if (M == 0) {
		cout << "Invalid input: Number of chains (M) should be greater than 0." << endl;
		return 1;
	}

	//generatedGraph = generateGraph(N, M);
	//Graph ReducedGraph = generateReducedGraph(chains, KeyNodes);


	
//#pragma omp parallel for
//	for (int i = 0; i < N; ++i) {
//		dijkstra(generatedGraph, i, N);
//	}
	

	//floydWarshall(generatedGraph.adjMatrix, N);
	//floydWarshall(ReducedGraph.adjMatrix, KeyNodes.size());



	//cout << "The adjacency matrix of the generated graph:" << endl;

	//读取数据后计算最短路
	Graph ReducedGraph(N);
	ReducedGraph.generateRandomWeights();

	if (writefile == 1) {

		writeAdjMatrixToCSV(to_string(N) + to_string(M) + "规模", ReducedGraph.adjMatrix);
	}
	if (readfile == 1) {
		ReducedGraph.adjMatrix = readAdjMatrixFromCSV(to_string(N) + to_string(M) + "规模");
	}
	//for (int i = 0; i < 1; ++i) {

	//}

	vector<vector<int>> ReducedNext=floydWarshall(ReducedGraph.adjMatrix, N);

	// 获取当前时间点
	auto start_time = std::chrono::high_resolution_clock::now();


	GraphResult ReducedResult(N,N, ReducedGraph.adjMatrix);

	cout << "Initializing get ODPaths!" << endl;

	//vector<vector<ODPaths>> ResultODPaths;
	std::vector<std::vector<ODPaths>> ResultODPaths(N, std::vector<ODPaths>(N, ODPaths(2)));

	cout << "Initializing all or nothing!" << endl;
	//1 初始化
	//1.1 初始化全有全无
	computeAllOrNothing(ResultODPaths, ReducedNext, N);

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


	//2 主循环
	int current_iteration = 0; // 当前迭代次数
	while (current_iteration < MaxIte) {
		std::cout << "Iteration: " << current_iteration << "-------------------------------------------------------------------------" << std::endl;
		//if (current_iteration % 100 == 0) {
		//	std::cout << "Iteration: " << current_iteration << std::endl;
		//}
		vector<vector<int>> loopReducedNext = floydWarshall(ReducedResult.bpr, N);

		std::cout << "Iteration: " << current_iteration <<" each OD "  <<std::endl;
//#pragma omp parallel for
		//2.1 对每个OD求解子问题
		for (int i = 0; i < N; ++i) {
			//if (i % 100 == 0) {
			//	std::cout << "Main loop i: " << i << std::endl;
			//}
			for (int j = 0; j < N; ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
				//if (j % 100 == 0) {
				//	std::cout << "Main loop j: " << j << std::endl;
				//}
				if (i != j) {
				 	//2.1.1 列生成
					//判断路径是否存在于 paths 中，不存在则添加到末尾
					//if (i == 99 && j == 33) {
					//	std::cout << "Iteration: " << current_iteration << " Adjustment " << i << "   " << j << "   " << std::endl;
					//}

					std::optional<Path> loopSP = getShortestPath(loopReducedNext, i, j);
					if (loopSP.has_value()) {
						//std::cout << "Optional contains a value." << std::endl;
						add_path_if_not_exists(ResultODPaths[i][j], loopSP.value());
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
						if (Adjustment> ResultODPaths[i][j].flow[path_id]){
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


		double total_convergence = calculate_total_convergence(ResultODPaths, ReducedResult.bpr);
		std::cout << "total_convergence: " << total_convergence << std::endl;
		if (total_convergence > 1 - ConvergenceAccuracy) { break; }
		// 增加当前迭代次数
		current_iteration++;



		// 可以在这里添加其他终止条件的判断逻辑
	}
	

	//


	// 获取当前时间点

	auto end_time = std::chrono::high_resolution_clock::now();

	// 计算时
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// 将时间转换为秒，并保留三位小数
	double durationInSeconds = duration.count() / 1000.0;

	// 输出执行时间
	std::cout << "程序执行时间: " << std::fixed << std::setprecision(3) << durationInSeconds << "秒" << std::endl;






	//printAdjacencyMatrix(generatedGraph);
	system("pause");
    return 0;
}

