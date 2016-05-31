#include <bits/stdc++.h>
#define MAXN 10
#define db(a) cout << #a << "=" << a << endl;
#define db2(a, b) cout << #a << "=" << a << " " << #b << " = " << b << endl;
#define db3(a, b, c) cout << #a << "=" << a << " " << #b << " = " << b << " " << #c << " = " << c << endl;
#define ROOT_PATH 0
#define NO_CITY -1
using namespace std;

int n, num_cities;
vector<vector<int> > G;

struct tour_t {
  int A[MAXN];
  int pos_city;
  int cost;
  tour_t(int ppos_city=0,int pcost=0){pos_city=ppos_city;cost=pcost;}
};
stack<int> s;

void mostrar_path (tour_t *tour);

tour_t *best_tour,*tour;
void Add_city(tour_t *tour,int next_city) {
  tour->A[tour->pos_city++] = next_city;
  if (tour->pos_city > 1) {
    tour->cost += G[tour->A[tour->pos_city - 2]][tour->A[tour->pos_city - 1]];
  }
  // db2(tour->city, w);
}

bool Best_tour(tour_t *tour) {
  bool res = best_tour->cost > tour->cost;
  return res;
}

bool Feasible(tour_t *tour, int next_city) {
  
  if (tour->pos_city == num_cities && next_city == ROOT_PATH) {
    return true;
  }
  for (int i = 0;i < tour->pos_city; i++) {
    if(tour->A[i] == next_city) return false;
  }
  
  return true;
}

void Remove_last_city(tour_t *tour) {
  if (tour->pos_city > 1) {
    tour->cost -= G[tour->A[tour->pos_city - 2]][tour->A[tour->pos_city - 1]];
    tour->A[tour->pos_city - 1] = 0;
    tour->pos_city--;
  }
}

void update_best_tour(tour_t *tour) {
  // memcpy(best_tour, tour, sizeof(tour_t));
  *best_tour = *tour;
}
int city_count(tour_t *tour) {
  return tour->pos_city;
}
void mostrar_path (tour_t *tour) {
  cout << "========= path ======" << endl;
  for (int i = 0; i < city_count(tour); i++) {
    cout << tour->A[i] << " ";
  }
  db(tour->cost)
  cout << endl;
}

int count_children(int city) {
    int c = 0;
    for (int i = 0; i < num_cities; i++) {
      if (G[city][i] != -1) c++;
    }
    return c;
}

int get_next_city(int cur_city, int nbr) {
  int next = 0;
  for (int i = 0; i < num_cities; i++) {
    if (G[cur_city][i] != -1) next++;
    if (next > nbr) return i;
  }
  return 0;
}

void TSP() {
  tour_t *tour = new tour_t;
  tour->cost = 0;
  tour->A[0] = ROOT_PATH;
  tour->pos_city = 1;
  
  s.push(ROOT_PATH);
  for (int city = 0; city < num_cities; city++) {
    if (ROOT_PATH == city) continue;
    s.push(city);
    while (s.size() != 1) {
      int cur_city = s.top();
      s.pop();
      if (cur_city == NO_CITY) {
        Remove_last_city(tour);
      } else {
        Add_city(tour, cur_city);
        if (city_count(tour) - 1 == num_cities) {
          if (Best_tour(tour)) {
            update_best_tour(tour);
            // mostrar_path(tour);
          }
          Remove_last_city(tour);
        } else {
          s.push(NO_CITY);
          for (int nbr = 0; nbr < count_children(cur_city); nbr++) {
            int next_city = get_next_city(cur_city, nbr);
            if (Feasible(tour, next_city)) {
              s.push(next_city);
            }
          }
        }
      }
    }
  }
  /* quitamos el primer elemento ROOT_PATH */
  s.pop();
}

void show_matrix () {
  for (int i = 0; i < num_cities; i++) {
    for (int j = 0; j < num_cities; j++) {
      cout << G[i][j] << "\t";
    }
    cout << endl;
  }
}
int main() {
  int times;
  scanf ("%d\n", &times);
  while (times--) {
    scanf("%d %d\n",&n, &num_cities);
    // db2(n, num_cities);
    // inicializar nodos
    G = vector<vector<int> >(num_cities, vector<int>(num_cities, -1));
    int a,b,w;
    // leyendo el grafo
    for(int i = 0;i < n; i++) {
      scanf("%d %d %d\n",&a,&b,&w);
      G[a][b] = w;
    }
    // show_matrix();
    best_tour = new tour_t;
    best_tour->cost = 1 << 30;
    TSP();
    printf("best_tour->cost %d\n", best_tour->cost);
    mostrar_path(best_tour);
  }
  return 0;
}