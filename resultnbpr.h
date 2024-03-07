#pragma once
// resultnbpr.h
#ifndef RESULTNBPR_H  // 防止头文件被多次引用
#define RESULTNBPR_H

#include <cmath>
#include <vector>
#include <chrono>
#include <utility> // for std::pair


const double big_m = 999999999999;



// 用户路径的结构
struct Path {
    std::vector<int> path;

    // 重载==运算符
    bool operator==(const Path& other) const {
        // 检查两个 vector 是否相等
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
    // 构造函数
    ODPaths(double d) : demand(d) {} // 使用成员初始化列表初始化 demand
};

// 函数定义：判断 Path 是否存在于 paths 中，不存在则添加到末尾
void add_path_if_not_exists(ODPaths& odPaths, const Path& path_to_add) {
    if (std::find(odPaths.paths.begin(), odPaths.paths.end(), path_to_add) == odPaths.paths.end()) {
        odPaths.paths.push_back(path_to_add);
        odPaths.flow.push_back(0);
    }
}


class GraphResult {
public:
    std::vector<std::vector<double>> flow; // 存储每条边上的流量
    std::vector<std::vector<double>> alpha; // 存储每条边上的 alpha
    std::vector<std::vector<double>> beta; // 存储每条边上的 beta
    std::vector<std::vector<int>> capacity; // 存储每条边上的能力
    std::vector<std::vector<double>> diff_bpr; // 存储每条边上的阻抗函数一阶导
    std::vector<std::vector<double>> bpr; // 存储每条边上的阻抗函数参数值
    std::vector<std::vector<int>> exceeded; // 存储每条边是否已达能力限制

    // 构造函数，根据点的数量创建对象
    GraphResult(int numPoints, int _capacity, const std::vector<std::vector<double>>& adjMatrix) {
        // 初始化 flow 向量大小为 numPoints
        flow.resize(numPoints, std::vector<double>(numPoints, 0));
        // 初始化 alpha 向量大小为 numPoints，并设置默认值为 0.15
        alpha.resize(numPoints, std::vector<double>(numPoints, 0.15));
        // 初始化 beta 向量大小为 numPoints，并设置默认值为 4
        beta.resize(numPoints, std::vector<double>(numPoints, 4));
        // 初始化 capacity 向量大小为 numPoints
        capacity.resize(numPoints, std::vector<int>(numPoints, _capacity));
        // 初始化 diff_bpr 向量大小为 numPoints，并设置默认值为 0
        diff_bpr.resize(numPoints, std::vector<double>(numPoints, 0));
        // 复制 adjMatrix 中的值给 bpr，并设置默认值为 距离
        bpr = adjMatrix;
        // 初始化 exceeded 向量大小为 numPoints，并设置默认值为 0
        exceeded.resize(numPoints, std::vector<int>(numPoints, 0));
    }

    // 其他方法和功能可以根据需要进行添加
};

double bpr_function(double t0, double alpha, double beta, double q, double qc) {
    return t0 * (1 + alpha * std::pow(q / qc, beta));
}

double bpr_function_derivative(double t0, double alpha, double beta, double q, double qc) {
    return t0 * alpha * beta * std::pow(q / qc, beta - 1) ;
}

// 计算某路径上所有路段 BPR 函数值的总和
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

// 计算某路径上所有路段 BPR 函数值的总和
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
            //cout << "i： " << path[i] << " j：" << path[i + 1] << endl;
            total_diff_bpr += diff_bpr[path[i]][path[i + 1]];
        }
    }
    return total_diff_bpr;
}



// 计算某OD路径与流量乘积的总和
double calculate_single_OD_total_flow_times_path(const std::vector<Path>& paths, const std::vector<double>& flow,
    const std::vector<std::vector<double>>& bpr) {
    double total = 0.0;
    for (int i = 0; i < paths.size(); ++i) {
        
        const Path& path = paths[i];
        double path_bpr = calculate_total_bpr(path.path,bpr); // 获得当前路径对应的bpr值
        // 计算当前路径流量乘以流量值的总和
        total += flow[i] * path_bpr; // 累加路径权重之和
        //cout << flow[i] << '\t' << path_bpr << '\t' << flow[i] * path_bpr << '\t';
    }
    return total;
}

// 计算单个 OD 所有路径对应的 BPR 函数值最小的一个
double calculate_min_bpr_for_od(const std::vector<Path>& paths,const std::vector<std::vector<double>>& bpr) {
    double min_bpr = std::numeric_limits<double>::max(); // 初始化为 double 类型的最大值

    for (int i = 0; i < paths.size(); ++i) {
        const Path& path = paths[i];
        double path_bpr = calculate_total_bpr(path.path, bpr); // 获得当前路径对应的bpr值
        // 如果当前路径的 BPR 函数值总和小于最小值，更新最小值
        if (path_bpr < min_bpr) {
            min_bpr = path_bpr;
        }
    }
    return min_bpr;
}



//获取path_bpr最小值及索引
std::pair<double, int> get_min_value_and_index(const std::vector<double>& values) {
    if (values.empty()) {
        // 如果向量为空，返回一个无效的最小值和索引
        return { std::numeric_limits<double>::quiet_NaN(), -1 };
    }

    double min_value = values[0]; // 初始化最小值为第一个元素
    int min_index = 0; // 初始化最小值的索引为0

    // 遍历向量，查找最小值及其索引
    for (int i = 1; i < values.size(); ++i) {
        if (values[i] < min_value) {
            min_value = values[i];
            min_index = i;
        }
    }

    return { min_value, min_index };
}

// 计算单个 OD 收敛性
double calculate_od_convergence(const std::vector<Path>& paths, const std::vector<double>& flow,
    const double demand, const std::vector<std::vector<double>>& bpr) {
    double min_bpr = demand * calculate_min_bpr_for_od(paths, bpr); // 计算收敛公式分子

    double total_flow_times_path = calculate_single_OD_total_flow_times_path(paths, flow, bpr); // 计算收敛公式分子分母

    // 返回收敛性指标（分子除以分母）
    return min_bpr / total_flow_times_path;
}

// 计算总的 ResultODPaths 收敛性
double calculate_total_convergence(const std::vector<std::vector<ODPaths>>& result_od_paths, const std::vector<std::vector<double>>& bpr) {
    double total_numerator = 0.0; // 总的分子
    double total_denominator = 0.0; // 总的分母

    // 遍历 ResultODPaths
#pragma omp parallel for reduction(+:total_numerator, total_denominator)
    for (int i = 0; i < result_od_paths.size(); ++i) {
        const auto& od_paths_vector = result_od_paths[i];
        for (int j = 0; j < od_paths_vector.size(); ++j) { // 从 i+1 开始遍历，跳过索引值相同的情况
            if (i != j) {
                const auto& od_paths = od_paths_vector[j];
                double min_bpr = calculate_min_bpr_for_od(od_paths.paths, bpr);
                double od_denominator = calculate_single_OD_total_flow_times_path(od_paths.paths, od_paths.flow, bpr);
                // 将当前 OD 的分子加到总的分子中
                total_numerator += od_paths.demand * min_bpr;
                //cout<< od_paths.demand<<'\t'<< min_bpr<<'\t'<< od_paths.demand * min_bpr<<endl;
                // 将当前 OD 的分母加到总的分母中
                total_denominator += od_denominator;
            }
        }
        if ((i+1) % 1000 == 0) {	std::cout  <<"i 分子： " << total_numerator << " 分母：" << total_denominator << endl;}
    }
    cout << "分子： "<< total_numerator<<" 分母："<< total_denominator << endl;
    // 返回总的收敛性指标
    return total_numerator / total_denominator;
}

//更新网络信息，flow、bpr、diff_bpr
void update_flow_bpr_derivative(const ODPaths& od_path,
    std::vector<std::vector<double>>& flow,
    std::vector<std::vector<double>>& diff_bpr,
    std::vector<std::vector<double>>& bpr,
    const std::vector<vector<double>> alpha,
    const std::vector<vector<double>> beta,
    const std::vector<vector<int>> capacity,
    const std::vector<vector<double>> adjMatrix) {
    // 获取当前时间点
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

    // 计算时
    auto _duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time - _start_time).count();


    // 输出执行时间
    std::cout << "函数程序1执行时间: " << _duration << "秒" << std::endl;




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

    // 计算时
    _duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time1 - _end_time).count();


    // 输出执行时间
    std::cout << "函数程序2执行时间: " <<  _duration << "秒" << std::endl;
     
}



#endif // RESULTNBPR_H