#pragma once
// resultnbpr.h
#ifndef RESULTNBPR_H  // ��ֹͷ�ļ����������
#define RESULTNBPR_H

#include <cmath>
#include <vector>
#include <chrono>
#include <utility> // for std::pair


const double big_m = 999999999999;



// �û�·���Ľṹ
struct Path {
    std::vector<int> path;

    // ����==�����
    bool operator==(const Path& other) const {
        // ������� vector �Ƿ����
        return std::equal(path.begin(), path.end(), other.path.begin(), other.path.end());
    }
};

template<typename T>
class Optional {
private:
    bool hasValue;
    T value;

public:
    Optional() : hasValue(false) {}
    Optional(const T& val) : hasValue(true), value(val) {}

    bool has_value() const { return hasValue; }
    const T& get_value() const {
        if (!hasValue) {
            throw std::runtime_error("Optional does not contain a value");
        }
        return value;
    }
};





struct ODPaths {
    std::vector<Path> paths;
    std::vector<double> flow;
    std::vector<double> path_bpr;
    std::vector<double> path_diff_bpr;

    const double demand;
    // ���캯��
    ODPaths(double d) : demand(d) {} // ʹ�ó�Ա��ʼ���б��ʼ�� demand
};

// �������壺�ж� Path �Ƿ������ paths �У�����������ӵ�ĩβ
void add_path_if_not_exists(ODPaths& odPaths, const Path& path_to_add) {
    if (std::find(odPaths.paths.begin(), odPaths.paths.end(), path_to_add) == odPaths.paths.end()) {
        odPaths.paths.push_back(path_to_add);
        odPaths.flow.push_back(0);
    }
}


class GraphResult {
public:
    std::vector<std::vector<double>> flow; // �洢ÿ�����ϵ�����
    std::vector<std::vector<double>> alpha; // �洢ÿ�����ϵ� alpha
    std::vector<std::vector<double>> beta; // �洢ÿ�����ϵ� beta
    std::vector<std::vector<int>> capacity; // �洢ÿ�����ϵ�����
    std::vector<std::vector<double>> diff_bpr; // �洢ÿ�����ϵ��迹����һ�׵�
    std::vector<std::vector<double>> bpr; // �洢ÿ�����ϵ��迹��������ֵ
    std::vector<std::vector<int>> exceeded; // �洢ÿ�����Ƿ��Ѵ���������

    // ���캯�������ݵ��������������
    GraphResult(int numPoints, int _capacity, const std::vector<std::vector<double>>& adjMatrix) {
        // ��ʼ�� flow ������СΪ numPoints
        flow.resize(numPoints, std::vector<double>(numPoints, 0));
        // ��ʼ�� alpha ������СΪ numPoints��������Ĭ��ֵΪ 0.15
        alpha.resize(numPoints, std::vector<double>(numPoints, 0.15));
        // ��ʼ�� beta ������СΪ numPoints��������Ĭ��ֵΪ 4
        beta.resize(numPoints, std::vector<double>(numPoints, 4));
        // ��ʼ�� capacity ������СΪ numPoints
        capacity.resize(numPoints, std::vector<int>(numPoints, _capacity));
        // ��ʼ�� diff_bpr ������СΪ numPoints��������Ĭ��ֵΪ 0
        diff_bpr.resize(numPoints, std::vector<double>(numPoints, 0));
        // ���� adjMatrix �е�ֵ�� bpr��������Ĭ��ֵΪ ����
        bpr = adjMatrix;
        // ��ʼ�� exceeded ������СΪ numPoints��������Ĭ��ֵΪ 0
        exceeded.resize(numPoints, std::vector<int>(numPoints, 0));
    }

    // ���������͹��ܿ��Ը�����Ҫ�������
};

double bpr_function(double t0, double alpha, double beta, double q, double qc) {
    return t0 * (1 + alpha * std::pow(q / qc, beta));
}

double bpr_function_derivative(double t0, double alpha, double beta, double q, double qc) {
    return t0 * alpha * beta * std::pow(q / qc, beta - 1) ;
}

// ����ĳ·��������·�� BPR ����ֵ���ܺ�
double calculate_total_bpr(const std::vector<int>& path, const std::vector<std::vector<double>>& bpr) {
    double total_bpr = 0.0;

    for (size_t i = 0; i < path.size() - 1; ++i) {
        int from = path[i];
        int to = path[i + 1];
        if (bpr[from][to] == numeric_limits < double > ::max()) { return big_m; }
        total_bpr += bpr[from][to];
    }
    return total_bpr;

    //double total_bpr_log = 0.0;

    //for (size_t i = 0; i < path.size() - 1; ++i) {
    //    int from = path[i];
    //    int to = path[i + 1];
    //    total_bpr_log += log(bpr[from][to]);
    //}

    //return exp(total_bpr_log);


}

// ����ĳ·��������·�� BPR ����ֵ���ܺ�
double calculate_total_diff_bpr(const std::vector<int>& path, const std::vector<int>& base_path, const std::vector<std::vector<double>>& diff_bpr) {
    double total_diff_bpr = 0.0;

    for (size_t i = 0; i < path.size() - 1; ++i) {
        int flag = 1;
        for (size_t j = 0; j < base_path.size() - 1; ++j) {
            if (path[i] == base_path[j] && path[i + 1] == base_path[j + 1]) {

                flag = 0;
                break;
            }
            //break;
        }
        if (flag == 1) {
            //cout << "i�� " << path[i] << " j��" << path[i + 1] << endl;
            total_diff_bpr += diff_bpr[path[i]][path[i + 1]];
        }
    }
    return total_diff_bpr;
}



// ����ĳOD·���������˻����ܺ�
double calculate_single_OD_total_flow_times_path(const std::vector<Path>& paths, const std::vector<double>& flow,
    const std::vector<std::vector<double>>& bpr) {
    double total = 0.0;
    for (int i = 0; i < paths.size(); ++i) {
        
        const Path& path = paths[i];
        double path_bpr = calculate_total_bpr(path.path,bpr); // ��õ�ǰ·����Ӧ��bprֵ
        // ���㵱ǰ·��������������ֵ���ܺ�
        total += flow[i] * path_bpr; // �ۼ�·��Ȩ��֮��
        //cout << flow[i] << '\t' << path_bpr << '\t' << flow[i] * path_bpr << '\t';
    }
    return total;
}

// ���㵥�� OD ����·����Ӧ�� BPR ����ֵ��С��һ��
double calculate_min_bpr_for_od(const std::vector<Path>& paths,const std::vector<std::vector<double>>& bpr) {
    double min_bpr = std::numeric_limits<double>::max(); // ��ʼ��Ϊ double ���͵����ֵ

    for (int i = 0; i < paths.size(); ++i) {
        const Path& path = paths[i];
        double path_bpr = calculate_total_bpr(path.path, bpr); // ��õ�ǰ·����Ӧ��bprֵ
        // �����ǰ·���� BPR ����ֵ�ܺ�С����Сֵ��������Сֵ
        if (path_bpr < min_bpr) {
            min_bpr = path_bpr;
        }
    }
    return min_bpr;
}



//��ȡpath_bpr��Сֵ������
std::pair<double, int> get_min_value_and_index(const std::vector<double>& values) {
    if (values.empty()) {
        // �������Ϊ�գ�����һ����Ч����Сֵ������
        return { std::numeric_limits<double>::quiet_NaN(), -1 };
    }

    double min_value = values[0]; // ��ʼ����СֵΪ��һ��Ԫ��
    int min_index = 0; // ��ʼ����Сֵ������Ϊ0

    // ����������������Сֵ��������
    for (int i = 1; i < values.size(); ++i) {
        if (values[i] < min_value) {
            min_value = values[i];
            min_index = i;
        }
    }

    return { min_value, min_index };
}

// ���㵥�� OD ������
double calculate_od_convergence(const std::vector<Path>& paths, const std::vector<double>& flow,
    const double demand, const std::vector<std::vector<double>>& bpr) {
    double min_bpr = demand * calculate_min_bpr_for_od(paths, bpr); // ����������ʽ����

    double total_flow_times_path = calculate_single_OD_total_flow_times_path(paths, flow, bpr); // ����������ʽ���ӷ�ĸ

    // ����������ָ�꣨���ӳ��Է�ĸ��
    return min_bpr / total_flow_times_path;
}

// �����ܵ� ResultODPaths ������
double calculate_total_convergence(const std::vector<std::vector<ODPaths>>& result_od_paths, const std::vector<std::vector<double>>& bpr) {
    double total_numerator = 0.0; // �ܵķ���
    double total_denominator = 0.0; // �ܵķ�ĸ

    // ���� ResultODPaths
#pragma omp parallel for reduction(+:total_numerator, total_denominator)
    for (int i = 0; i < result_od_paths.size(); ++i) {
        const auto& od_paths_vector = result_od_paths[i];
        for (int j = 0; j < od_paths_vector.size(); ++j) { // �� i+1 ��ʼ��������������ֵ��ͬ�����
            if (i != j) {
                const auto& od_paths = od_paths_vector[j];
                double min_bpr = calculate_min_bpr_for_od(od_paths.paths, bpr);
                double od_denominator = calculate_single_OD_total_flow_times_path(od_paths.paths, od_paths.flow, bpr);
                // ����ǰ OD �ķ��Ӽӵ��ܵķ�����
                total_numerator += od_paths.demand * min_bpr;
                //cout<< od_paths.demand<<'\t'<< min_bpr<<'\t'<< od_paths.demand * min_bpr<<endl;
                // ����ǰ OD �ķ�ĸ�ӵ��ܵķ�ĸ��
                total_denominator += od_denominator;
            }
        }
        if ((i+1) % 1000 == 0) {	std::cout  <<"i ���ӣ� " << total_numerator << " ��ĸ��" << total_denominator << endl;}
    }
    cout << "���ӣ� "<< total_numerator<<" ��ĸ��"<< total_denominator << endl;
    // �����ܵ�������ָ��
    return total_numerator / total_denominator;
}

//����������Ϣ��flow��bpr��diff_bpr
void update_flow_bpr_derivative(const ODPaths& od_path,
    std::vector<std::vector<double>>& flow,
    std::vector<std::vector<double>>& diff_bpr,
    std::vector<std::vector<double>>& bpr,
    const std::vector<vector<double>> alpha,
    const std::vector<vector<double>> beta,
    const std::vector<vector<int>> capacity,
    const std::vector<vector<double>> adjMatrix) {
    // ��ȡ��ǰʱ���
    auto _start_time = std::chrono::high_resolution_clock::now();

    for (int k = 0; k < od_path.paths.size(); ++k) {
        const Path& path = od_path.paths[k];
        double singel_flow = od_path.flow[k];
        for (int t = 0; t < path.path.size() - 1; ++t) {
            int from = path.path[t];
            int to = path.path[t + 1];
            flow[from][to] += singel_flow;
        }
    }



    auto _end_time = std::chrono::high_resolution_clock::now();

    // ����ʱ
    auto _duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time - _start_time).count();


    // ���ִ��ʱ��
    std::cout << "��������1ִ��ʱ��: " << _duration << "��" << std::endl;




    for (int k = 0; k < od_path.paths.size(); ++k) {
        const Path& path = od_path.paths[k];
        double singel_flow = od_path.flow[k];
        for (int t = 0; t < path.path.size() - 1; ++t) {
            int from = path.path[t];
            int to = path.path[t + 1];
            bpr[from][to] = bpr_function(adjMatrix[from][to], alpha[from][to], beta[from][to], flow[from][to], capacity[from][to]);
            diff_bpr[from][to] = bpr_function(adjMatrix[from][to], alpha[from][to], beta[from][to], flow[from][to], capacity[from][to]);
        }
    }


    auto _end_time1 = std::chrono::high_resolution_clock::now();

    // ����ʱ
    _duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time1 - _end_time).count();


    // ���ִ��ʱ��
    std::cout << "��������2ִ��ʱ��: " <<  _duration << "��" << std::endl;
     
}



#endif // RESULTNBPR_H