#include <bits/stdc++.h>
#define MAXN 100
#define db(a) cout << #a << "=" << a << endl;
#define db2(a, b) cout << #a << "=" << a << " " << #b << " = " << b << endl;
#define ROOT_PATH 0
using namespace std;
vector<pair<int,int> > G[MAXN];
int n, num_cities;

struct tour_t {
  int A[MAXN];
  int city;
  int cost;
  tour_t(int pcity=0,int pcost=0){city=pcity;cost=pcost;}
};

tour_t *best_tour,*tour;
void Add_city(tour_t *tour,int city, int w) {
  tour->A[tour->city++] = city;
  tour->cost += w;
}

bool Best_tour(tour_t *tour) {
  return best_tour->cost > tour->cost;
}

bool Feasible(tour_t *tour,int city) {
  for(int i=0;i < tour->city;i++) {
    if (tour->A[i] == city) return false;
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
void TSP(tour_t *tour) {
  if(city_count(tour) == num_cities) {
    int w = 0, idx = -1,   = tour->A[tour->city - 1];
    for (int i = 0; i < G[cur_city].size(); i++) {
      if (G[cur_city][i].first == ROOT_PATH) {
        idx = i; break;
      }
    }
    w = G[cur_city][idx].second;
    Add_city(tour, ROOT_PATH, w);
    if(Best_tour(tour)) { 
      update_best_tour(tour);
    }
    Remove_last_city(tour, ROOT_PATH, w);
  } else {
    int city = tour->A[tour->city - 1];
    for(int nbr = 0; nbr < G[city].size(); nbr++) {
      int cur_city = G[city][nbr].first;
      int w = G[city][nbr].second;
      if(Feasible(tour, cur_city)) {
      	Add_city(tour, cur_city, w);
      	TSP(tour);
      	Remove_last_city(tour,cur_city, w);
      }
    }
  }
}

int main() {
  scanf("%d %d\n",&n, &num_cities);
  // inicializar nodos
  fill(G,G+n,vector<pair<int,int> >());
  int a,b,w;
  // leyendo el grafo
  for(int i=0;i<n;i++) {
    scanf("%d %d %d\n",&a,&b,&w);
    G[a].push_back(make_pair(b,w));
  }
  
  tour = new tour_t;
  best_tour = new tour_t;
  best_tour->cost = 1 << 30;
  tour->cost = 0;
  tour->A[0] = ROOT_PATH;
  tour->city = 1;
  
  TSP(tour);
  printf("best_tour->cost %d\n", best_tour->cost);
  
  return 0;
}