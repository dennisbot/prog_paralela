#include <bits/stdc++.h>
#include <pthread.h>
#define MAXN 100
#define db(a) cout << #a << "=" << a << endl;
#define db2(a, b) cout << #a << "=" << a << " " << #b << " = " << b << endl;
#define dbs(a, out) pthread_mutex_lock(&out); db(a); pthread_mutex_unlock(&out);
#define ROOT_PATH 0
#define thread_count 4
using namespace std;

pthread_t *thread_handles;
pthread_mutex_t lock, left_stack, term_mutex, out;
pthread_cond_t term_cond_var;

vector<pair<int,int> > G[MAXN];
int n, num_cities, threads_in_cond_wait, my_stack_size;
struct tour_t {
  int A[MAXN];
  int city;
  int cost;
  tour_t(int pcity=0,int pcost=0){city=pcity;cost=pcost;}
};

vector<tour_t*> subtours;
stack<tour_t*> new_stack;

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

stack<tour_t*> Split_my_stack(stack<tour_t*> my_stack) {
  int my_stack_size = my_stack.size();
  stack<tour_t*> new_stack;
  for (int i = 0; i < my_stack_size / 2; i++) {
    new_stack.push(my_stack.top());
    my_stack.pop();
  }
  return new_stack;
}

bool Terminated(stack<tour_t *> &my_stack) {
  my_stack_size = my_stack.size();
  if (my_stack_size >= 2 && threads_in_cond_wait > 0 && new_stack.size() == 0) {
    pthread_mutex_lock(&term_mutex);
    if (threads_in_cond_wait > 0 && new_stack.size() == 0) {
      new_stack = Split_my_stack(my_stack); 
      pthread_cond_signal(&term_cond_var);
      pthread_mutex_unlock(&term_mutex);
      return 0; /* Terminated = false; don’t quit */ 
    } 
    else if (!my_stack.empty()) { /* Keep working */
        return 0; /* Terminated = false; don’t quit */ 
    } 
    else { /* My_stack is empty */
      pthread_mutex_unlock(&term_mutex);
      if (threads_in_cond_wait == thread_count - 1) {
        /* Last thread running */
        threads_in_cond_wait++;
        pthread_cond_broadcast(&term_cond_var); 
        pthread_mutex_unlock(&term_mutex);
        return 1; /* Terminated = true; quit */
      }  
      else { /* Other threads still working, wait for work */ 
        threads_in_cond_wait++;
        while (pthread_cond_wait(&term_cond_var, &term_mutex) != 0); /* We’ve been awakened */
        if (threads_in_cond_wait < thread_count) { /* We got work */ 
          my_stack = new_stack;
          stack<tour_t*> ntour;
          swap(new_stack, ntour);
          threads_in_cond_wait--;
          pthread_mutex_unlock(&term_mutex);
          return 0; 
        }/* Terminated = false */  
        else { /* All threads done */
          pthread_mutex_unlock(&term_mutex);
          return 1; /* Terminated = true; quit */ 
        }
      } /* else wait for work */ 
    } /* else my_stack is empty */
  }
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
// void *TSP_DFS(tour_t* tour) {
void *TSP_DFS(void* param_p) {
  tour_t* my_tour = (tour_t*) param_p;
  pthread_mutex_lock(&out);
  printf("rank thread = %d of a Total of %d Threads With PThreads Static\n", -1, thread_count);
  mostrar_path(my_tour);
  pthread_mutex_unlock(&out);
  // return NULL;

  stack<tour_t*> l_stack;
  l_stack.push(my_tour);
  
  if (!new_stack.empty()) {
    pthread_mutex_lock(&left_stack);
      while (!new_stack.empty()) {
        l_stack.push(new_stack.top());
        new_stack.pop();
      }
    pthread_mutex_unlock(&left_stack);
  }
  
  while (!Terminated(l_stack)) {
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
        pthread_mutex_lock(&lock);
        if(Best_tour(cur_tour)) {
          // mostrar_path(cur_tour);
          update_best_tour(cur_tour);
        }
        pthread_mutex_unlock(&lock);
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
  long thread;
  scanf ("%d\n", &times);
  while (times--) {
    subtours.clear();
    stack<tour_t*> ntour;
    swap(new_stack, ntour);

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
    
    while (!q.empty()) {
      subtours.push_back(q.front().first);
      q.pop();
    }
    for (int i = 0; i < subtours.size(); i++) {
      mostrar_path(subtours[i]);
    }
    
    threads_in_cond_wait = 0;
    thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
    pthread_mutex_init(&lock, NULL);
    pthread_mutex_init(&left_stack, NULL);
    pthread_mutex_init(&out, NULL);
    pthread_mutex_init(&term_mutex, NULL);
    pthread_cond_init(&term_cond_var, NULL);
    
    db(subtours.size())
    db(thread_count);
    
    for (int i = thread_count; i < subtours.size(); i++) {
      new_stack.push(subtours[i]);
    }

    for (thread = 0; thread < thread_count && thread < subtours.size(); thread++) {
      // pthread_mutex_lock(&out);
      // db(thread)
      // pthread_mutex_unlock(&out);
      pthread_create(&thread_handles[thread], NULL, TSP_DFS, (void*) subtours[thread]);
    }
    for (thread = 0; thread < thread_count && thread < subtours.size(); thread++) {
      pthread_join(thread_handles[thread], NULL);
    }

    pthread_mutex_destroy(&lock);
    pthread_mutex_destroy(&left_stack);
    pthread_mutex_destroy(&out);
    pthread_cond_destroy(&term_cond_var);

    // el arreglo que usas es qq
    printf("best_tour->cost %d\n", best_tour->cost);
  }
  return 0;
}