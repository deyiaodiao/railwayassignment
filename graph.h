#pragma once
// graph.h
#ifndef GRAPH_H  // 防止头文件被多次引用
#define GRAPH_H

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <limits>

#include <random>
#include <stack>
#include <unordered_map>
#include <string>
#include <fstream>
#include <algorithm> // 包含 std::find

using namespace std;


//const double INF =99999999; // 无穷大值


class Graph {
public:
	vector<vector<double>> adjMatrix; // 存储权重

public:
	Graph(int numNodes) {
		adjMatrix.resize(numNodes, vector<double>(numNodes, INF));
	}

	void addEdge(int u, int v, double weight) {
		adjMatrix[u][v] = weight;
		adjMatrix[v][u] = weight; // For undirected graph
	}
	// 获取权重
	double getWeight(int u, int v) {
		return adjMatrix[u][v];
	}
	void printGraph() {
		for (int i = 0; i < adjMatrix.size(); ++i) {
			std::cout << "Node " << i << " connected to: ";
			for (int j = 0; j < adjMatrix[i].size(); ++j) {
				if (adjMatrix[i][j] != 0) {
                    std::cout << j << " (weight: " << adjMatrix[i][j] << ") ";
				}
			}
            std::cout << endl;
		}
	}

    // Dijkstra算法
    double shortestPath(int start, int end, vector<int>& path) {
        int numNodes = adjMatrix.size();
        if (start >= numNodes || end >= numNodes) {
            std::cout << "Error: Node index out of range." << endl;
            return -1; // 返回错误值
        }
        vector<double> dist(numNodes, INF);
        vector<bool> visited(numNodes, false);
        vector<int> pred(numNodes, -1); // 前驱节点数组，用于记录路径
        dist[start] = 0;

        for (int i = 0; i < numNodes - 1; ++i) {
            int minIndex = -1;
            double minDist = INF;
            for (int j = 0; j < numNodes; ++j) {
                if (!visited[j] && dist[j] < minDist) {
                    minDist = dist[j];
                    minIndex = j;
                }
            }

            if (minIndex == -1 || minIndex == end)
                break;

            visited[minIndex] = true;

            for (int j = 0; j < numNodes; ++j) {
                if (!visited[j] && adjMatrix[minIndex][j] != INF && dist[minIndex] + adjMatrix[minIndex][j] < dist[j]) {
                    dist[j] = dist[minIndex] + adjMatrix[minIndex][j];
                    pred[j] = minIndex; // 更新前驱节点数组
                }
            }
        }

        // 从终点回溯到起点，构建路径
        stack<int> s;
        int current = end;
        while (current != start) {
            s.push(current);
            current = pred[current];
        }
        s.push(start);

        // 将路径信息存储到 path 数组中
        while (!s.empty()) {
            path.push_back(s.top());
            s.pop();
        }

        return dist[end];
    }


    void generateRandomWeights() {
        // 使用随机数生成器引擎
        default_random_engine generator(time(nullptr));
        // 定义一个分布，生成浮点数在 [1, 100] 范围内
        uniform_real_distribution<double> distribution(1.0, 100.0);
#pragma omp parallel for
        for (int i = 0; i < adjMatrix.size(); ++i) {
            for (int j = i + 1; j < adjMatrix[i].size(); ++j) { // 只需要设置上三角
                double weight = distribution(generator);
                adjMatrix[i][j] = weight;
                adjMatrix[j][i] = weight; // 由于是无向图，对称设置权重
            }
            adjMatrix[i][i] = 0;
        }
    }

};

// 边的信息结构体，用于存储节点与距离
struct EdgeInfo {
	int node_id;
	double distance;
};

// 节点类
class Node {
public:
	// 构造函数
	Node(int identifier) :
		identifier(identifier){}

	// 成员变量
	int identifier;                             // 节点的唯一标识符
	int chain_number=-1;                           // 节点所处链的编号

	int identifier_B=-1;                             // 节点的唯一标识符
	std::unordered_map<int, EdgeInfo> edges_A;  // 在网络A中与其他节点的连接信息
	std::unordered_map<int, EdgeInfo> edges_B;  // 在网络B中与其他节点的连接信息
};






#endif // GRAPH_H

//// 新建一个节点，假设唯一标识符为1001，所处链的编号为1
//Node newNode(1001, 1);
//
//// 修改节点在网络A中的连接信息
//newNode.edges_A[1002] = { 1002, 10.5f };  // 与节点1002连接，距离为10.5
//newNode.edges_A[1003] = { 1003, 20.0f };  // 与节点1003连接，距离为20.0
//
//// 修改节点在网络B中的连接信息
//newNode.edges_B[2001] = { 2001, 15.3f };  // 与节点2001连接，距离为15.3
//newNode.edges_B[2002] = { 2002, 25.8f };  // 与节点2002连接，距离为25.8
//
//// 将节点添加到全局变量的vector中
//nodes.push_back(newNode);
//
//// 输出节点的信息
//std::cout << "Node Info:\n";
//std::cout << "Identifier: " << newNode.identifier << "\n";
//std::cout << "Chain Number: " << newNode.chain_number << "\n";
//
//// 输出节点在网络A中的连接信息
//std::cout << "Edges in Network A:\n";
//for (const auto& edge : newNode.edges_A) {
//	std::cout << "  To Node: " << edge.second.node_id << ", Distance: " << edge.second.distance << "\n";
//}
//
//// 输出节点在网络B中的连接信息
//std::cout << "Edges in Network B:\n";
//for (const auto& edge : newNode.edges_B) {
//	std::cout << "  To Node: " << edge.second.node_id << ", Distance: " << edge.second.distance << "\n";
//}
//
