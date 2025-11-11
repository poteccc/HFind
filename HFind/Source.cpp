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
    int max_edge = 0;
    for (int i = 0; i < Graph.size(); i++) if (Graph[i][Graph[i].size() - 1] > max_edge) max_edge = Graph[i][Graph[i].size() - 1];
    for (int i = 0; i <= max_edge; i++)Graph_Transposed.push_back(vector<int>());
    for (int i = 0; i < Graph.size(); i++) for (int j = 0; j < Graph[i].size(); j++) Graph_Transposed[Graph[i][j]].push_back(i);
    for (int i = 0; i <= Graph_Transposed.size(); i++) sort(Graph_Transposed.begin(), Graph_Transposed.end());
    return Graph_Transposed;
}
vector<int> g(vector<vector<int>> Graph, vector<int> input, bool is_vertex = true) { // работает правильно
    int massize = input.size();
    vector<int> ans;
    vector<int> b;
    for (int i = 0; i < input.size(); i++)b.push_back(0);
    if (!is_vertex) {
        while (true) {
            bool flag = 1;
            for (int j = 0; j < input.size(); j++) {
                if (b[j] == Graph[input[j]].size()) return ans;
                if (Graph[input[j]][b[j]] > Graph[input[0]][b[0]]) {
                    b[0]++;
                    flag = 0;
                }
                if (Graph[input[j]][b[j]] < Graph[input[0]][b[0]]) {
                    b[j]++;
                    flag = 0;
                }
            }
            if (flag) {
                ans.push_back(Graph[input[0]][b[0]]);
                b[0]++;
            }
        }
    }
    else {
        auto GraphT = Graph_Transpose(Graph);
        while (true) {
            bool flag = 1;
            for (int j = 0; j < input.size(); j++) {
                if (b[j] == GraphT[input[j]].size()) return ans;
                if (GraphT[input[j]][b[j]] > GraphT[input[0]][b[0]]) {
                    b[0]++;
                    flag = 0;
                }
                if (GraphT[input[j]][b[j]] < GraphT[input[0]][b[0]]) {
                    b[j]++;
                    flag = 0;
                }
            }
            if (flag) {
                ans.push_back(GraphT[input[0]][b[0]]);
                b[0]++;
            }
        }
    }
    return ans;
}
vector<int> gg(vector<vector<int>> Graph, vector<int> input, bool is_vertex = true) { // работает правильно
    if (is_vertex) return g(Graph, g(Graph, input, true), false);
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
            delta_level[Xu].push_back(i);
            if (v[level_num - 1] == delta_U[i] - 1) p--;
            else p = level_num - 1;
            if (p >= 0) for (int j = level_num - 1; j >= p; j--) v[j] = v[p] + j - p + 1;
        }
    }
}
void HFindMCS(map<vector<int>, vector<int>>& MCS, vector<vector<int>> Graph) { // работает правильно
    int delta = 0, * delta_U = new int[Graph.size()];
    for (int i = 0; i < Graph.size(); i++) {
        delta_U[i] = Graph[i].size();
        if (delta < Graph[i].size()) delta = Graph[i].size();
    }
    map<vector<int>, vector<int>>* levels = new map<vector<int>, vector<int>>[delta];
    for (int i = delta; i > 0; i--) {
        Clu(levels[i - 1], Graph, i, delta_U);
        for (auto it = levels[i - 1].begin(); it != levels[i - 1].end(); ++it) {
            bool flag = 1;
            if (!MCS.count(it->second)) MCS[it->second] = it->first;
            //D-вершины, C-рёбра
           /*cout << endl << "C:";           //Вывод MCS
           for (int i = 0; i < it->first.size(); i++) {
            cout << it->first[i] << ' ';
           }
           cout << endl << "D:";
           for (int i = 0; i < it->second.size(); i++) {
            cout <<
            it->second[i] << ' ';
    }
    cout << endl;
    //*/
        }
    }
    /*int SIZE = 0;                //Вывод числа комбинаций
    for (int i = 0; i < delta; i++) SIZE += levels[i].size();
    cout << SIZE;*/
}
/*void Dynamic_MCS(vector<pair<vector<int>, vector<int>>>& MCS, vector<int> new_line, vector<vector<int>> old_graph) {
 //vector<pair<vector<int>, vector<int>>> new_vertex_combinations;
 vector<vector<int>> new_vertex_combinations;
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
   for (int j = 0; j < level_num; j++) Xu.push_back(new_line[level_num]);
   new_vertex_combinations.push_back(Xu);
   if()
   if (v[level_num - 1] == delta_U[level_num] - 1) p--;
   else p = level_num - 1;
   if (p >= 0) for (int j = level_num - 1; j >= p; j--) v[j] = v[p] + j - p + 1;
  }
 }
 /*
 HFindMCS(new_vertex_combinations, graph);
 for (int i = 0; i < new_vertex_combinations.size(); i++) {
  bool flag = 0;
  vector<int> C = gg(old_Graph, new_vertex_combinations[i].first, true);      //хз что это, не помню
  for (int j = 0; j < MCS.size(); j++) {
   if (MCS[j].first == C) {
    MCS[j].second.push_back(old_Graph.size());
    flag = 1;
    break;
   }
  }
  if (!flag) {
   vector<int>D = g(old_Graph, new_vertex_combinations[i].first, true);
   MCS.push_back(pair<vector<int>, vector<int>>(C, D));
  }
 }*/
 //}
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
    /*vector<vector<int>> graph(6), graph2;
    graph[0].push_back(0);
    graph[1].push_back(0);
    graph[1].push_back(1);
    graph[2].push_back(0);
    graph[2].push_back(1);
    graph[2].push_back(2);
    graph[3].push_back(3);
    graph[4].push_back(2);
    graph[5].push_back(2);
    graph2 = graph;
    Graph_output(Graph_Transpose(graph));*/
    vector<vector<int>> graph;
    const int N = 10, M = 10, delta = 5;
    //cout << N / delta * 2 + 1;
    for (int i = 0; i < N; i++) {                                //рандомная генерация графа
        graph.push_back(vector<int>());
        int k = 0;
        for (int j = 0; j < M; j++) {
            if (rand() % (N / delta + 1) == 0 && k < delta) {
                graph[i].push_back(j);
                k++;
            }
        }
    }
    Graph_output(graph);
    map<vector<int>, vector<int>> MCS;
    HFindMCS(MCS, graph);
    //cout << "____________________________________________________________________________________________";
    cout << MCS.size();
    //Dynamic_MCS(MCS, a, graph);
    for (int i = 0; i < MCS.size(); i++) {
        cout << endl << "C:";
        for (int j = 0; j < MCS[i].first.size(); j++) {     //Вывод MCS
            cout << MCS[i].first[j] << ' ';
        }
        cout << endl << "D:";
        for (int j = 0; j < MCS[i].second.size(); j++) {
            cout << MCS[i].second[j] << ' ';
        }
        cout << endl;
    }//*/
}