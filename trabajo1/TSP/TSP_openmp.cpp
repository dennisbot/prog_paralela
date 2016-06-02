#include <bits/stdc++.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#define MAXN 100
#define db(a) cout << #a << "=" << a << endl;
#define db2(a, b) cout << #a << "=" << a << " " << #b << " = " << b << endl;
#define ROOT_PATH 0
#define thread_count 4
using namespace std;
vector<pair<int,int> > G[MAXN];
int n, num_cities;
struct tour_t {
  int A[MAXN];
  int city;
  int cost;
  tour_t(int pcity=0,int pcost=0){city=pcity;cost=pcost;}
};
queue<pair<tour_t*, int> > q;
tour_t *best_tour,*tour;
void Add_city(tour_t *tour,int city, int w) {
  tour->A[tour->city++] = city;
  tour->cost += w;
  // db2(tour->city, w);
}

bool Best_tour(tour_t *tour) {
  bool res = best_tour->cost > tour->cost;
  return res;
}

bool Feasible(tour_t *tour, int city) {
  for(int i=0;i < tour->city;i++) {
    // db2(tour->A[i], city);
    if(tour->A[i] == city) return false;
  }
  return true;
}

void Remove_last_city(tour_t *tour,int city, int w) {
  tour->A[tour->city - 1] = 0;
  tour->cost -= w;
  tour->city--;
}

void update_best_tour(tour_t *tour) {
  memcpy(best_tour, tour, sizeof(tour_t));
}
int city_count(tour_t *tour) {
  return tour->city;
}
void mostrar_path (tour_t *tour) {
  cout << "========= best ======" << endl;
  for (int i = 0; i < city_count(tour); i++) {
    cout << tour->A[i] << " ";
  }
  cout << endl;
  cout << tour->cost << endl;

}
queue<pair<tour_t*, int> > TSP_BFS(int nivel) {

  tour_t *tour = new tour_t;
  tour->cost = 0;
  tour->A[0] = ROOT_PATH;
  tour->city = 1;
  
  q.push(make_pair(tour, 0));
  int times = 0;
  while (!q.empty()) {
    int q_size = q.size();
    pair<tour_t*, int> ptour_t = q.front();
    tour_t *cur_tour = ptour_t.first;
    int depth = ptour_t.second;
    if (depth >= nivel) break;
    // NUM_THREADS
    q.pop();
    // seguimos expandiendo en BFS
    int city = cur_tour->A[cur_tour->city - 1];
    for (int nbr = 0; nbr < G[city].size(); nbr++) {
      int cur_city = G[city][nbr].first;
      int w = G[city][nbr].second;
      if(Feasible(cur_tour, cur_city)) {
        Add_city(cur_tour, cur_city, w);
        tour_t *new_tour = new tour_t;
        memcpy(new_tour, cur_tour, sizeof(tour_t));
        q.push(make_pair(new_tour, depth + 1));
        Remove_last_city(cur_tour, cur_city, w);
      }
    }
  }
  return q;
}
void TSP_DFS(tour_t* tour) {
  printf("rank thread = %d of a Total of %d Threads With OpenMP\n", omp_get_thread_num(), thread_count);
  
  stack<tour_t*> l_stack;
  l_stack.push(tour);

  while (!l_stack.empty()) {
    tour_t *cur_tour = l_stack.top();
    l_stack.pop();
    if (city_count(cur_tour) == num_cities) {
      int w = 0, idx = -1, cur_city = cur_tour->A[cur_tour->city - 1];
      for (int i = 0; i < G[cur_city].size(); i++) {
        if (G[cur_city][i].first == ROOT_PATH) {
          idx = i; break;
        }
      }
      w = G[cur_city][idx].second;
      Add_city(cur_tour, ROOT_PATH, w);
      
      if(Best_tour(cur_tour)) {
        #pragma omp critical
        {
          if(Best_tour(cur_tour)) {
            // mostrar_path(cur_tour);
            update_best_tour(cur_tour);
          }
        }
      }  
      Remove_last_city(cur_tour, ROOT_PATH, w);
    }
    else {
      int city = cur_tour->A[cur_tour->city - 1];
      
      for (int nbr = 0; nbr < G[city].size(); nbr++) {
        int cur_city = G[city][nbr].first;
        int w = G[city][nbr].second;
        if(Feasible(cur_tour, cur_city)) {
          Add_city(cur_tour, cur_city, w);
          tour_t *new_tour = new tour_t;
          memcpy(new_tour, cur_tour, sizeof(tour_t));
          l_stack.push(new_tour);
          Remove_last_city(cur_tour,cur_city, w);
        }
      }
    }
  }
}

int main() {
  int times;
  scanf ("%d\n", &times);
  while (times--) {
    scanf("%d %d\n",&n, &num_cities);
    // db2(n, num_cities);
    // inicializar nodos
    fill(G, G + n, vector<pair<int,int> >());
    int a, b, w;
    // leyendo el grafo
    for(int i = 0; i < n; i++) {
      scanf("%d %d %d\n", &a, &b, &w);
      G[a].push_back(make_pair(b, w));
    }
    best_tour = new tour_t;
    best_tour->cost = 1 << 30;
    int nivel = num_cities / 2;
    queue<pair<tour_t*, int> > q = TSP_BFS(nivel);
    vector<tour_t*> subtours;
    while (!q.empty()) {
      subtours.push_back(q.front().first);
      q.pop();
    }
    for (int i = 0; i < subtours.size(); i++) {
      mostrar_path(subtours[i]);
    }
    #pragma omp parallel for num_threads(thread_count)
    for (int i = 0; i < subtours.size(); i++) {
      TSP_DFS(subtours[i]);
    }

    // el arreglo que usas es qq
    printf("best_tour->cost %d\n", best_tour->cost);
  }
  return 0;
}