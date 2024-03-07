#pragma once
// graph.h
#ifndef GRAPH_H  // ��ֹͷ�ļ����������
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
#include <algorithm> // ���� std::find

using namespace std;


//const double INF =99999999; // �����ֵ


class Graph {
public:
	vector<vector<double>> adjMatrix; // �洢Ȩ��

public:
	Graph(int numNodes) {
		adjMatrix.resize(numNodes, vector<double>(numNodes, INF));
	}

	void addEdge(int u, int v, double weight) {
		adjMatrix[u][v] = weight;
		adjMatrix[v][u] = weight; // For undirected graph
	}
	// ��ȡȨ��
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

    // Dijkstra�㷨
    double shortestPath(int start, int end, vector<int>& path) {
        int numNodes = adjMatrix.size();
        if (start >= numNodes || end >= numNodes) {
            std::cout << "Error: Node index out of range." << endl;
            return -1; // ���ش���ֵ
        }
        vector<double> dist(numNodes, INF);
        vector<bool> visited(numNodes, false);
        vector<int> pred(numNodes, -1); // ǰ���ڵ����飬���ڼ�¼·��
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
                    pred[j] = minIndex; // ����ǰ���ڵ�����
                }
            }
        }

        // ���յ���ݵ���㣬����·��
        stack<int> s;
        int current = end;
        while (current != start) {
            s.push(current);
            current = pred[current];
        }
        s.push(start);

        // ��·����Ϣ�洢�� path ������
        while (!s.empty()) {
            path.push_back(s.top());
            s.pop();
        }

        return dist[end];
    }


    void generateRandomWeights() {
        // ʹ�����������������
        default_random_engine generator(time(nullptr));
        // ����һ���ֲ������ɸ������� [1, 100] ��Χ��
        uniform_real_distribution<double> distribution(1.0, 100.0);
#pragma omp parallel for
        for (int i = 0; i < adjMatrix.size(); ++i) {
            for (int j = i + 1; j < adjMatrix[i].size(); ++j) { // ֻ��Ҫ����������
                double weight = distribution(generator);
                adjMatrix[i][j] = weight;
                adjMatrix[j][i] = weight; // ����������ͼ���Գ�����Ȩ��
            }
            adjMatrix[i][i] = 0;
        }
    }

};

// �ߵ���Ϣ�ṹ�壬���ڴ洢�ڵ������
struct EdgeInfo {
	int node_id;
	double distance;
};

// �ڵ���
class Node {
public:
	// ���캯��
	Node(int identifier) :
		identifier(identifier){}

	// ��Ա����
	int identifier;                             // �ڵ��Ψһ��ʶ��
	int chain_number=-1;                           // �ڵ��������ı��

	int identifier_B=-1;                             // �ڵ��Ψһ��ʶ��
	std::unordered_map<int, EdgeInfo> edges_A;  // ������A���������ڵ��������Ϣ
	std::unordered_map<int, EdgeInfo> edges_B;  // ������B���������ڵ��������Ϣ
};






#endif // GRAPH_H

//// �½�һ���ڵ㣬����Ψһ��ʶ��Ϊ1001���������ı��Ϊ1
//Node newNode(1001, 1);
//
//// �޸Ľڵ�������A�е�������Ϣ
//newNode.edges_A[1002] = { 1002, 10.5f };  // ��ڵ�1002���ӣ�����Ϊ10.5
//newNode.edges_A[1003] = { 1003, 20.0f };  // ��ڵ�1003���ӣ�����Ϊ20.0
//
//// �޸Ľڵ�������B�е�������Ϣ
//newNode.edges_B[2001] = { 2001, 15.3f };  // ��ڵ�2001���ӣ�����Ϊ15.3
//newNode.edges_B[2002] = { 2002, 25.8f };  // ��ڵ�2002���ӣ�����Ϊ25.8
//
//// ���ڵ���ӵ�ȫ�ֱ�����vector��
//nodes.push_back(newNode);
//
//// ����ڵ����Ϣ
//std::cout << "Node Info:\n";
//std::cout << "Identifier: " << newNode.identifier << "\n";
//std::cout << "Chain Number: " << newNode.chain_number << "\n";
//
//// ����ڵ�������A�е�������Ϣ
//std::cout << "Edges in Network A:\n";
//for (const auto& edge : newNode.edges_A) {
//	std::cout << "  To Node: " << edge.second.node_id << ", Distance: " << edge.second.distance << "\n";
//}
//
//// ����ڵ�������B�е�������Ϣ
//std::cout << "Edges in Network B:\n";
//for (const auto& edge : newNode.edges_B) {
//	std::cout << "  To Node: " << edge.second.node_id << ", Distance: " << edge.second.distance << "\n";
//}
//
