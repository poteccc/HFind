#include <iostream>
#include <math.h>
#include <ctime>
#include <array>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <omp.h>
#include <list>
#include <Windows.h>
#include <algorithm>
using namespace std;
const int num_threads = 4;
template<typename T>
int bin_search(vector<T> mas, T search) { // работает правильно
    int left = 0, right = mas.size(), step = 0;
    while (right - left > 1) {
        step = (right + left) / 2;
        if (mas[step] <= search) left = step;
        if (mas[step] >= search) right = step;
    }
    if (mas[step] == search) return step;
    else return -1;
}

vector<vector<int>> Graph_Transpose(vector<vector<int>> Graph) { // работает правильно
    vector<vector<int>> Graph_Transposed;
    //for (int i = 0; i < M; i++)Graph_Transposed.push_back(vector<int>());
    //int* mas_b = new int[N];
    int max_edge = 0;
    for (int i = 0; i < Graph.size(); i++) if (Graph[i][Graph[i].size() - 1] > max_edge) max_edge = Graph[i][Graph[i].size() - 1];
    for (int i = 0; i <= max_edge; i++)Graph_Transposed.push_back(vector<int>());
    for (int i = 0; i < Graph.size(); i++) for (int j = 0; j < Graph[i].size(); j++) Graph_Transposed[Graph[i][j]].push_back(i);
    return Graph_Transposed;
}

vector<int> g(vector<vector<int>> Graph, vector<int> input, bool is_edges = true) { // работает правильно
    vector<int> ans;
    if (is_edges) {
        for (int i = 0; i < Graph.size(); i++) {
            if (Graph[i].size() < input.size()) continue;
            int j = 0, k = 0;
            while (j < Graph[i].size() && k < input.size()) {
                if (Graph[i][j] == input[k]) {
                    j++;
                    k++;
                    continue;
                }
                if (Graph[i][j] < input[k]) {
                    j++;
                    continue;
                }
                if (Graph[i][j] > input[k]) {
                    break;
                }
            }
            if (k == input.size()) ans.push_back(i);
        }
    }
    else {
        auto GraphT = Graph_Transpose(Graph);
        for (int i = 0; i < GraphT.size(); i++) {
            if (GraphT[i].size() < input.size()) continue;
            int j = 0, k = 0;
            while (j < GraphT[i].size() && k < input.size()) {
                if (GraphT[i][j] == input[k]) {
                    j++;
                    k++;
                    continue;
                }
                if (GraphT[i][j] < input[k]) {
                    j++;
                    continue;
                }
                if (GraphT[i][j] > input[k]) {
                    break;
                }
            }
            if (k == input.size()) ans.push_back(i);
        }
    }
    return ans;
}

vector<int> gg(vector<vector<int>> Graph, vector<int> input, bool is_edges = true) { // работает правильно
    if (is_edges) return g(Graph, g(Graph, input, true), false);
    else return g(Graph, g(Graph, input, false), true);
}

void Clu(map<vector<int>, vector<int>>& delta_level, vector<vector<int>> Graph, int level_num, int* delta_U) { // работает правильно
    for (int i = 0; i < Graph.size(); i++) {
        if (level_num > delta_U[i]) continue;
        if (level_num == delta_U[i]) {
            delta_level[Graph[i]].push_back(i);
            continue;
        }
        int p = level_num - 1;
        vector<int> v;
        for (int t = 0; t < level_num; t++) v.push_back(t);
        while (p >= 0) {
            vector<int> Xu;
            for (int j = 0; j < level_num; j++) Xu.push_back(Graph[i][v[j]]);
            delta_level[Xu].push_back(i);                                       //Ключ - рёбра(номера столбцов), значение - вершины(номера строк)
            if (v[level_num - 1] == delta_U[i] - 1) p--;
            else p = level_num - 1;
            if (p >= 0) for (int j = level_num - 1; j >= p; j--) v[j] = v[p] + j - p + 1;
        }
    }
}

void Clu_parallel(map<vector<int>, vector<int>>& delta_level, vector<vector<int>> Graph, int level_num, int* delta_U) {
    map<vector<int>, vector<int>> delta_thread[num_threads];
    
#pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < Graph.size(); i++) {
        if (level_num > delta_U[i]) continue;
        if (level_num == delta_U[i]) {
            delta_thread[omp_get_thread_num()][Graph[i]].push_back(i);
            continue;
        }
        int p = level_num - 1;
        vector<int> v;
        for (int t = 0; t < level_num; t++) v.push_back(t);
        while (p >= 0) {
            vector<int> Xu;
            for (int j = 0; j < level_num; j++) Xu.push_back(Graph[i][v[j]]);
            delta_thread[omp_get_thread_num()][Xu].push_back(i);
            if (v[level_num - 1] == delta_U[i] - 1) p--;
            else p = level_num - 1;
            if (p >= 0) for (int j = level_num - 1; j >= p; j--) v[j] = v[p] + j - p + 1;
        }
    }
    for (int p = 0; p < num_threads; p++) {
        for (auto iter = delta_thread[omp_get_thread_num()].begin(); iter != delta_thread[omp_get_thread_num()].end(); iter++) {
            if (delta_level.count(iter->first)) {
                //P[k - 1][Xu] = P[k - 1][Xu] + char(U[i] + 1);
                for (auto itr = iter->first.begin(); itr != iter->first.end(); itr++) delta_level[iter->first].push_back(*itr);
            }
            else {
                delta_level[iter->first] = iter->second;
            }
        }
        delta_thread[p].clear();
    }
}

void HFindMCS(vector<pair<vector<int>, vector<int>>>& MCS, vector<vector<int>> Graph,int parallel = 0) { // работает правильно
    int delta = 0, * delta_U = new int[Graph.size()];
    for (int i = 0; i < Graph.size(); i++) {
        delta_U[i] = Graph[i].size();
        if (delta < Graph[i].size()) delta = Graph[i].size();
    }
    map<vector<int>, vector<int>>* levels = new map<vector<int>, vector<int>>[delta];
    for (int i = delta; i > 0; i--) {
        if (parallel == 0) Clu(levels[i - 1], Graph, i, delta_U);
        if (parallel == 1) Clu_parallel(levels[i - 1], Graph, i, delta_U);
        for (auto it = levels[i - 1].begin(); it != levels[i - 1].end(); ++it) {
            bool flag = 1;
            for (int j = 0; j < MCS.size();j++) {
                if (MCS[j].second == it->second) {
                    flag = 0;
                    break;
                }
            }
            if (flag) MCS.push_back(make_pair(it->first, it->second));
        }
    }
}

void HFindMCS_parallel(vector<pair<vector<int>, vector<int>>>& MCS, vector<vector<int>> Graph) {
    int delta = 0, * delta_U = new int[Graph.size()];
    for (int i = 0; i < Graph.size(); i++) {
        delta_U[i] = Graph[i].size();
        if (delta < Graph[i].size()) delta = Graph[i].size();
    }
    map<vector<int>, vector<int>>* levels = new map<vector<int>, vector<int>>[delta];
#pragma omp parallel for num_threads(num_threads)
    for (int i = delta; i > 0; i--) Clu(levels[i - 1], Graph, i, delta_U);
#pragma omp barrier
    for (int i = delta; i > 0; i--) {
        //Clu(levels[i - 1], Graph, i, delta_U);
        for (auto it = levels[i - 1].begin(); it != levels[i - 1].end(); ++it) {
            bool flag = 1;
            for (int j = 0; j < MCS.size(); j++) {
                if (MCS[j].second == it->second) {
                    flag = 0;
                    break;
                }
            }
            if (flag) MCS.push_back(make_pair(it->first, it->second));
        }
    }
    /*int SIZE = 0;                //Вывод числа комбинаций
    for (int i = 0; i < delta; i++) SIZE += levels[i].size();
    cout << SIZE;*/
}

void Dynamic_MCS(vector<pair<vector<int>, vector<int>>>& MCS, vector<int> new_line, vector<vector<int>> old_graph) {
    //vector<pair<vector<int>, vector<int>>> new_vertex_combinations;
    vector<vector<int>> new_vertex_combinations;
    vector<int> new_vertex_combination, new_edges_combination;
    vector<vector<int>> graph;
    graph.push_back(new_line);
    int* delta_U = new int;
    *delta_U = new_line.size();
    for (int level_num = 0; level_num < new_line.size(); level_num++) {
        int p = level_num - 1;
        vector<int> v;
        for (int t = 0; t < level_num; t++) v.push_back(t);
        while (p >= 0) {
            vector<int> Xu;
            for (int j = 0; j < level_num; j++) Xu.push_back(new_line[v[j]]);
            new_vertex_combination = g(old_graph, Xu);
            bool flag = 1;
            for (int i = 0; i < new_vertex_combinations.size(); i++) {
                if (new_vertex_combinations[i] == new_vertex_combination) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                new_edges_combination = g(old_graph, new_vertex_combination, false);
                bool flag = 1;
                for (int j = 0; j < MCS.size(); j++) {
                    if (MCS[j].first == new_edges_combination) {
                        flag = 0;
                        break;
                    }
                }
                if (flag) {
                    MCS.push_back(make_pair(new_edges_combination, new_vertex_combination));
                }
            }
            if (v[level_num - 1] == delta_U[level_num] - 1) p--;
            else p = level_num - 1;
            if (p >= 0) for (int j = level_num - 1; j >= p; j--) v[j] = v[p] + j - p + 1;
        }
    }
}

void Graph_output(vector<vector<int>> Graph) { // работает правильно
    int max_edge = 0;
    for (int i = 0; i < Graph.size(); i++) if (Graph[i][Graph[i].size() - 1] > max_edge) max_edge = Graph[i][Graph[i].size() - 1];
    for (int i = 0; i < Graph.size(); i++) {
        int k = 0;
        for (int j = 0; j < Graph[i].size(); j++) {
            while (k++ < Graph[i][j]) cout << 0 << " ";
            cout << 1 << " ";
        }
        while (k++ <= max_edge) cout << 0 << " ";
        cout << endl;
    }
}

int main() {
    srand(time(0));
    vector<vector<int>> graph, new_graph;
    const int N = 5, M = 5, delta = 5;

    for (int i = 0; i < N; i++) {                                //рандомная генерация графа
        graph.push_back(vector<int>());
        int k = 0;
        for (int j = 0; j < M; j++) {
            if (rand() % (N / delta + 1) == 0 && k < delta) {
                graph[i].push_back(j);
                k++;
            }
        }
        if (k == 0) graph[i].push_back(rand() % N);
    }

    new_graph = graph;
    vector<int> new_vertex;
    int k = 0;

    for (int i = 0; i < M; i++) {
        if (rand() % (N / delta + 1) == 0 && k < delta) {
            new_vertex.push_back(i);
            k++;
        }
    }

    new_graph.push_back(new_vertex);
    cout << "na4alo" << endl;
    Graph_output(graph);
    cout << endl;
    int k = 0;
    for (int j = 0; j < new_vertex.size(); j++) {
        while (k++ < new_vertex[j]) cout << 0 << " ";
        cout << 1 << " ";
    }
    while (k++ <= M) cout << 0 << " ";

    //Graph_output(Graph_Transpose(graph));
    vector<pair<vector<int>, vector<int>>> MCS1;
    vector<pair<vector<int>, vector<int>>> MCS2;
    int time = omp_get_wtime();
    HFindMCS(MCS1, graph);
    Dynamic_MCS(MCS1, new_vertex, graph);
    cout << endl << omp_get_wtime() - time << endl << MCS1.size();
    time = omp_get_wtime();
    //HFindMCS_parallel(MCS2, graph);
    HFindMCS(MCS2, new_graph);
    cout << endl << omp_get_wtime() - time << endl << MCS1.size();
}