#define INT_MAX 1000
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include "fiboheap.h"
#include "argparse.hpp"
#include <string>
#include <queue>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <fstream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

typedef std::pair<double, int> dij_Pair;

int minDistance(std::vector<double> &dist, std::vector<bool> &sptSet, long int V_m)
{
   // Initialize min value
   double min = INT_MAX, min_index;

   for (int v = 0; v < V_m; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;

   return min_index;
}

double triangle_area_calculater(double a, double b, double c)
{
  double p = (a + b + c)/2;
  return sqrt(p*(p-a)*(p-b)*(p-c));
}

double l2_norm(Eigen::MatrixXd &V, long int i, long int j)
{
  return sqrt(pow(V(i,0) - V(j,0), 2) + pow(V(i,1) - V(j,1), 2) +
    pow(V(i,2) - V(j,2), 2));
}

struct dij_out
{
  double min_dist;
  std::vector<int> prev;
  std::vector<double> dist;
  std::vector<int> parent;
};



dij_out dijkstra(std::vector<std::vector<dij_Pair>> adj_vec, long int src, long int des)
{


    long int V_m = adj_vec.size();
    // std::cout << "V_m is " << V_m << std::endl;
    std::vector<double> dist(V_m);     // The output array.  dist[i] will hold the shortest
    std::vector<int> prev(V_m);
    std::vector<bool> sptSet(V_m); // sptSet[i] will be true if vertex i is included in shortest
    // distance from src to i
    std::vector<int> parent;

    // path tree or shortest distance from src to i is finalized

    // Initialize all distances as INFINITE and stpSet[] as false
    for (size_t i = 0; i < V_m; i++)
    dist[i] = INT_MAX, sptSet[i] = false, prev[i] = -1;

    // Distance of source vertex from itself is always 0
    dist[src] = 0;

    // Find shortest path for all vertices
    for (int count = 0; count < V_m-1; count++)
    {
      // std::cout << "count number " << count << std::endl;
      // Pick the minimum distance vertex from the set of vertices not
      // yet processed. u is always equal to src in the first iteration.
      int u = minDistance(dist, sptSet, V_m);
      // if(u == des)
      //   break;
      // std::cout << "min vertex is " << u << std::endl;

      // Mark the picked vertex as processed
      sptSet[u] = true;

      // Update dist value of the adjacent vertices of the picked vertex.
      std::vector<dij_Pair>::iterator i;
      for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++){

      // Update dist[v] only if is not in sptSet, there is an edge from
      // u to v, and total weight of path from src to  v through u is
      // smaller than current value of dist[v]
      long int v = (*i).second;
      int weight = (*i).first;
      if (!sptSet[v] && dist[u] != INT_MAX
      && dist[u]+weight < dist[v])
      {
        dist[v] = dist[u] + weight;
        prev[v] = u;
      }
      }
    }


    int k = des;
    if(prev[k]>=0 || k==src)
    {
      while(k>=0)
      {
        parent.push_back(k);
        k = prev[k];
      }
    }

    dij_out retval;
    retval.min_dist = dist[des];
    retval.prev = prev;
    retval.dist = dist;
    retval.parent = parent;

    return retval;
}


dij_out dijkstra_min_heap(std::vector<std::vector<dij_Pair>> adj_vec, int src, int des)
{
      int V_m = adj_vec.size();
      // std::cout << "V_m is " << V_m << std::endl;
      std::vector<double> dist(V_m);     // The output array.  dist[i] will hold the shortest
      std::vector<int> prev(V_m);
      std::vector<bool> sptSet(V_m); // sptSet[i] will be true if vertex i is included in shortest
                      // distance from src to i
     std::vector<int> parent;
     std::priority_queue <dij_Pair, std::vector<dij_Pair>, std::greater<dij_Pair> > pq;

                     // path tree or shortest distance from src to i is finalized
     // Distance of source vertex from itself is always 0
     dist[src] = 0;
     pq.push(std::make_pair(0, src));

     // Initialize all distances as INFINITE and stpSet[] as false
     for (size_t i = 0; i < V_m; i++)
     {
       if(i != src)
       {
         pq.push(std::make_pair(INT_MAX, i)), dist[i] = INT_MAX;
       }
       prev[i] = -1;
     }

     // Find shortest path for all vertices
     while (!pq.empty())
     {
       // std::cout << "count number " << count << std::endl;
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = pq.top().second;
       pq.pop();
       // std::cout << "min vertex is " << u << std::endl;
       // if(u == des)
       //   break;

       std::vector<dij_Pair>::iterator i;
       // Update dist value of the adjacent vertices of the picked vertex.
       for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++)
       {
         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
        long int v = (*i).second;
         double weight = (*i).first;

         // If there is shorted path to v through u.
         if (dist[v] > dist[u] + weight)
         {
             // Updating distance of v
             dist[v] = dist[u] + weight;
             pq.push(std::make_pair(dist[v], v));
             prev[v] = u;
         }
       }
     }


     int u = des;
     if(prev[u]>=0 || u==src)
     {
       while(u>=0)
       {
         parent.push_back(u);
         u = prev[u];
       }
     }

     dij_out retval;
     retval.min_dist = dist[des];
     retval.prev = prev;
     retval.dist = dist;
     retval.parent = parent;

     return retval;

     // print the constructed distance array
     // printSolution(dist, V_m);
}

dij_out dijkstra_fibo_heap(std::vector<std::vector<dij_Pair>> adj_vec, int src, int des)
{
      int V_m = adj_vec.size();
      // std::cout << "V_m is " << V_m << std::endl;
     std::vector<double> dist(V_m);     // The output array.  dist[i] will hold the shortest
     std::vector<int> prev(V_m);
     std::vector<bool> sptSet(V_m); // sptSet[i] will be true if vertex i is included in shortest

                      // distance from src to i
     std::vector<int> parent;
     FibHeap <dij_Pair> pq;

                     // path tree or shortest distance from src to i is finalized
     // Distance of source vertex from itself is always 0

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V_m; i++)
     {
       if(i == src)
       {
         pq.push(std::make_pair(0, i));
         dist[i] = 0;
       }
       else
       {
         pq.push(std::make_pair(INT_MAX, i)), dist[i] = INT_MAX;
       }
       prev[i] = -1;
     }

     // Find shortest path for all vertices
     while (!pq.empty())
     {
       // std::cout << "count number " << count << std::endl;
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = pq.top().second;
       pq.pop();
       // std::cout << "min vertex is " << u << std::endl;
       // if(u == des)
       //   break;

       std::vector<dij_Pair>::iterator i;
       // Update dist value of the adjacent vertices of the picked vertex.
       for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++)
       {
         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         int v = (*i).second;
         double weight = (*i).first;

         // If there is shorted path to v through u.
         if (dist[v] > dist[u] + weight)
         {
             // Updating distance of v
             dist[v] = dist[u] + weight;
             pq.push(std::make_pair(dist[v], v));
             prev[v] = u;
         }
       }
     }

     int u = des;
     if(prev[u]>=0 || u==src)
     {
       while(u>=0)
       {
         parent.push_back(u);
         u = prev[u];
       }
     }

     dij_out retval;
     retval.min_dist = dist[des];
     retval.prev = prev;
     retval.dist = dist;
     retval.parent = parent;

     return retval;
}

std::vector<double> iso_curve_signature(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
   std::vector<std::vector<dij_Pair>> &my_pair, std::vector<int> &fps, int index)
{
  // start the geodesic iso-curve signature
  std::cout << "Start the iso function" << '\n';
    long int iso_start = fps[index];
    int k_iso = 20;
    std::vector<double> iso_curve_bin(k_iso, 0);
    dij_out temp_dij = dijkstra_fibo_heap(my_pair, iso_start, 0);
    double max_geodesic = *max_element(temp_dij.dist.begin(), temp_dij.dist.end());
    double k_unit = max_geodesic/k_iso;
    for (size_t i = 0; i < F.rows(); i++) {
      long int a = F(i,0), b = F(i,1), c = F(i,2);
      double k_a = temp_dij.dist[a];
      double k_b = temp_dij.dist[b];
      double k_c = temp_dij.dist[c];
      int radii_number = round((k_a + k_b + k_c)/(3*k_unit));
      double r_i = k_unit*radii_number;
      // std::cout << "The iso-curve bin number is " << radii_number << '\n';
      if((k_a < r_i && k_b > r_i && k_c > r_i) || (k_a > r_i && k_b < r_i && k_c < r_i))
      {
        double a1, a2;
        Eigen::VectorXd p1, p2;
        a1 = abs(r_i - k_a)/abs(k_b - k_a);
        a2 = abs(r_i - k_a)/abs(k_c - k_a);
        p1 = (1 - a1)*V.row(a) + a1*V.row(b);
        p2 = (1 - a2)*V.row(a) + a2*V.row(c);
        iso_curve_bin[radii_number] += (p1 - p2).norm();
      }
      else if((k_b < r_i && k_a > r_i && k_c > r_i) || (k_b > r_i && k_a < r_i && k_c < r_i))
      {
        double a1, a2;
        Eigen::VectorXd p1, p2;
        a1 = abs(r_i - k_b)/abs(k_a - k_b);
        a2 = abs(r_i - k_b)/abs(k_c - k_b);
        p1 = (1 - a1)*V.row(b) + a1*V.row(a);
        p2 = (1 - a2)*V.row(b) + a2*V.row(c);
        iso_curve_bin[radii_number] += (p1 - p2).norm();
      }
      else if((k_c < r_i && k_a > r_i && k_b > r_i) || (k_c > r_i && k_a < r_i && k_b < r_i))
      {
        double a1, a2;
        Eigen::VectorXd p1, p2;
        a1 = abs(r_i - k_c)/abs(k_a - k_c);
        a2 = abs(r_i - k_c)/abs(k_b - k_c);
        p1 = (1 - a1)*V.row(c) + a1*V.row(a);
        p2 = (1 - a2)*V.row(c) + a2*V.row(b);
        iso_curve_bin[radii_number] += (p1 - p2).norm();
      }
    }
    return iso_curve_bin;
}

std::vector<double> bilateral_descriptor(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
   std::vector<std::vector<dij_Pair>> &my_pair, std::vector<int> &fps,
   long int index, long int fin)
{
  // Initialize the bilateral_hist with zeros
  int k_iso = 20;
  std::vector<double> bilateral_hist(k_iso, 0);
  // std::cout << "Min path distance is " << my_out.min_dist << '\n';
  // std::vector<int> roi_elementes;
  long int iso_start = fps[index];
  std::vector<bool> is_roi_element(V.rows(), false);
  size_t count_roi = 0;
  dij_out my_out = dijkstra_fibo_heap(my_pair, iso_start, fin);
  std::vector<int> parent_m = my_out.parent;
  double p_q_distance = my_out.min_dist;
  double distance_threshold = 0.8*my_out.min_dist;
  double bilateral_unit = p_q_distance/20;
  Eigen::MatrixXd path_geodesic(parent_m.size(), V.rows());
  for (size_t i = 0; i < parent_m.size(); i++) {
    dij_out my_temp_out = dijkstra_fibo_heap(my_pair, parent_m[parent_m.size()-i-1], fin);
    for (size_t j = 0; j < V.rows(); j++) {
      path_geodesic(i,j) = my_temp_out.dist[j];
    }
  }
  for (size_t i = 0; i < V.rows(); i++) {
    double min_geo_to_path = INT_MAX;
    for (size_t j = 0; j < parent_m.size(); j++) {
      int k = parent_m[j];
      if(path_geodesic(j, i) < min_geo_to_path && i!=k)
      {
        min_geo_to_path = path_geodesic(j, i);
      }
    }
    if(min_geo_to_path < distance_threshold &&
      path_geodesic(0, i) < p_q_distance && path_geodesic(parent_m.size()-1, i) < p_q_distance)
    {
      // roi_elementes.push_back(i);
      is_roi_element[i] = true;
      count_roi++;
      // int k_class = round(min_geo_to_path/bilateral_unit);
    }
  }

  // std::cout << "ROI has " << count_roi << " elements" << '\n';
  // iterate over the faces and compute bilateral histogram
  for (size_t i = 0; i < F.rows(); i++) {
    int a = F(i,0), b = F(i,1), c = F(i,2);
    if (is_roi_element[a] && is_roi_element[b] && is_roi_element[c]) {
      double middle_bi = (path_geodesic(0, a) + path_geodesic(0, b)
        + path_geodesic(0, c))/3;
      int k_class = round(middle_bi/bilateral_unit);
      double a_b = l2_norm(V, a, b);
      double b_c = l2_norm(V, b, c);
      double a_c = l2_norm(V, a, c);
      double face_area = triangle_area_calculater(a_b, a_b, a_b);
      bilateral_hist[k_class] += face_area;
      // std::cout << "Z value is " << Z(i) << '\n';
    }
  }
  return bilateral_hist;
}


int main(int argc, char *argv[])
{
  // ArgumentParser parser;
  // parser.addArgument("-f", "--file");
  // parser.addArgument("-s", "--source", 0, true);
  // parser.addArgument("-d", "--destination", 50, true);

  // parse the command-line arguments - throws if invalid format


  // Load a mesh in OFF format
  igl::readOFF(argv[1], V, F);

  std::string start, finish, fps_opt;
  start = argv[2];
  finish = argv[3];
  int st = 0;
  int fin = 0;
  int N_samples = 10;
  long int m_v = V.rows();
  int n_v = V.cols();
  long int m_f = F.rows();
  int n_f = F.cols();
  bool calc_n_n = false;

  // Eigen::MatrixXd V_cost(m_v, m_v);
  // Eigen::MatrixXi V_adj(m_v, m_v);
  // V_cost = Eigen::MatrixXd::Zero(m_v, m_v);
  // V_adj = Eigen::MatrixXi::Zero(m_v, m_v);
  // std::vector<std::vector<int> > VF;
  // std::vector<std::vector<int> > VFi;

  try
  {
    st = std::stoi(start);
    fin = std::stoi(finish);
    N_samples = std::stoi(argv[4]);
  }
  catch(...)
  {
    std::cout << "Input vertices cannot be converted to int" <<
    std::endl;
  }

  if(st >= m_v)
  {
    std::cerr << "Source vertex index cannot be larger than the vertex number: " << m_v-1 << '\n';
    std::exit(1);
  }
  else if(fin > m_v )
  {
    std::cerr << "Destination vertex index cannot be larger than the vertex number: " << m_v-1 << '\n';
    std::exit(1);
  }


  std::cout << "The start is " << st << " and finish is " <<
    fin << std::endl;
  std::cout << "The V size is " << m_v << '\n';
  std::cout << "The F size is " << m_f << '\n';
  // for(size_t i=0; i<m_v; i++)
  // {
  //   for(size_t j=0; j<m_v; j++)
  //   {
  //     double length = sqrt(pow(V(i,0) - V(j,0), 2) + pow(V(i,1) - V(j,1), 2) +
  //       pow(V(i,2) - V(j,2), 2));
  //       V_cost(i,j) = length;
  //   }
  // }

  std::vector<std::vector<double>> A;
  igl::adjacency_list(F,A, true);
  std::vector<std::vector<dij_Pair>> my_pair(A.size());
  std::cout << "Adjacency list computed" << '\n';

    // for(size_t j=0; j<A[0].size(); j++)
    // {
    //   std::cout << A[st][j] << " ";
    // }
    // std::cout << std::endl;

  for(size_t i=0; i<A.size(); i++)
  {
    for(size_t j=0; j<A[i].size(); j++)
    {
      long int k = A[i][j];
      double weight = l2_norm(V, i, k);
      // V_cost(i, k) = weight;
      my_pair[i].push_back(std::make_pair(weight, k));
    }
  }

  std::cout << "my_pair computed" << '\n';

  // for(size_t i=0; i<10; i++)
  // {
  //   // std::cout << "The results are \n";
  //   for(size_t j=0; j<10; j++)
  //   {
  //     std::cout << V_cost(i,j) << " ";
  //   }
  //   std::cout << std::endl;
  // }

  double min_dist;
  double max_geodesic = 0;
  dij_out my_out;

  if(calc_n_n)
  {
    std::clock_t start_time = std::clock();
    double duration;
    Eigen::MatrixXd V_geodesic(m_v, m_v);

    for (size_t i = 0; i < m_v; i++) {
      std::cout << "Dijkstra step " << i << '\n';
      dij_out my_out_temp = dijkstra_fibo_heap(my_pair, i, fin);
      if(i == st)
        my_out = my_out_temp;
      for (size_t j = 0; j < m_v; j++) {
        // std::cout << "Dijkstra second step " << j << '\n';
        V_geodesic(i,j) = my_out_temp.dist[j];
        if(my_out_temp.dist[j] > max_geodesic && i == st)
          max_geodesic = my_out_temp.dist[j];
      }
    }

    duration = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;
    std::cout << "Operation took "<< duration << "seconds for array dij" << std::endl;
    // my_out = dijkstra(my_pair, st, fin);
    // my_out = dijkstra_min_heap(my_pair, st, fin);
    // my_out = dijkstra_fibo_heap(my_pair, st, fin);

    std::ofstream myfile_n_n;
    myfile_n_n.open ("first_geodesic.txt");
    for (size_t i = 0; i < m_v; i++) {
      for (size_t j = 0; j < m_v; j++) {
        myfile_n_n << V_geodesic(i,j) << " ";
      }
      myfile_n_n << "\n";
    }
    myfile_n_n.close();
  }
  else
  {
    std::clock_t start_time = std::clock();
    double duration;

    my_out = dijkstra(my_pair, st, fin);

    duration = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;
    std::cout << "Operation took "<< duration << "seconds for array dij" << std::endl;

    start_time = std::clock();

    my_out = dijkstra_min_heap(my_pair, st, fin);

    duration = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;
    std::cout << "Operation took "<< duration << "seconds for min heap dij" << std::endl;

    start_time = std::clock();

    my_out = dijkstra_fibo_heap(my_pair, st, fin);

    duration = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;
    std::cout << "Operation took "<< duration << "seconds for fibo heap dij" << std::endl;

    for (size_t i = 0; i < my_out.dist.size(); i++) {
      if(my_out.dist[i] > max_geodesic)
        max_geodesic = my_out.dist[i];
    }
  }


  std::cout << "The minimum distance is " << my_out.min_dist << std::endl;
  std::cout << "The maximum distance is " << max_geodesic << std::endl;

  int u = fin;
  std::vector<int> parent_m = my_out.parent;
  // if(my_out.prev[u]>=0 || u==st)
  // {
  //   while(u>=0)
  //   {
  //     parent_m.push_back(u);
  //     u = my_out.prev[u];
  //   }
  // }
// print the route nodes
  for(size_t i=0; i<parent_m.size(); i++)
  {
	  std::cout << parent_m[i] << " ";
 }
 std::cout << std::endl;

 Eigen::MatrixXd path_geodesic(parent_m.size(), m_v);
 std::vector<double> max_geodesic_vec(parent_m.size(), 0);
 for (size_t i = 0; i < parent_m.size(); i++) {
   dij_out temp_path = dijkstra_fibo_heap(my_pair, parent_m[parent_m.size()-i-1], fin);
   for (size_t j = 0; j < m_v; j++) {
     path_geodesic(i,j) = temp_path.dist[j];
     if(temp_path.dist[j] > max_geodesic_vec[i])
      max_geodesic_vec[i] = temp_path.dist[j];
   }
 }

// put the coordinates of the vertices in matrices to draw the line
Eigen::MatrixXd P1(parent_m.size()-1, 3);
Eigen::MatrixXd P2(parent_m.size()-1, 3);
for(size_t i=0; i<parent_m.size()-1; i++)
{
  int vertex1 = parent_m[i];
  int vertex2 = parent_m[i+1];
  double epsilon = 0.0001;
  // adding 0.1 to coordinates to make the line visible
  P1(i,0) = V(vertex1,0) + epsilon, P1(i,1) = V(vertex1,1) + epsilon, P1(i,2) = V(vertex1,2) + epsilon;
  P2(i,0) = V(vertex2,0) + epsilon, P2(i,1) = V(vertex2,1) + epsilon, P2(i,2) = V(vertex2,2) + epsilon;
}

// this the FPS part

bool isVertexUsed[m_v];
Eigen::MatrixXd fps_geo_dist(N_samples, m_v);
for (size_t i = 0; i < m_v; i++) {
  isVertexUsed[i] = false;
}
isVertexUsed[0] = true;
std::vector<int> fps(N_samples);
fps[0] = st;
for (size_t i = 1; i < N_samples; i++) {
  // int len_try = m_v-i;
  double max = 0;
  int max_place = 0;
  dij_out my_dij = dijkstra_fibo_heap(my_pair, fps[i-1], 0);
  for (size_t j = 0; j < m_v; j++) {
    fps_geo_dist(i-1, j) = my_dij.dist[j];
  }
  // std::vector<double> distances(len_try);
  for (size_t j = 0; j < m_v; j++) {
    if(isVertexUsed[j])
      continue;

    double min_inner = INT_MAX;
    for (size_t k = 0; k < i; k++) {
      // dij_out my_dij = dijkstra(V_cost, fps[k], j);
      size_t vert_index = fps[k];
      double curr_dist = fps_geo_dist(k, j);
      if(curr_dist < min_inner)
        min_inner = curr_dist;
    }
    if(max < min_inner)
    {
      max = min_inner;
      max_place = j;
    }

  }
  fps[i] = max_place;
  isVertexUsed[max_place] = true;
  std::cout << "The iteration number is " << i << '\n';
  std::cout << "The selected vertex is " << max_place << '\n';
}
// put the fps points into to P for displaying them
Eigen::MatrixXd P(N_samples, 3);
for(size_t i=0; i<N_samples; i++)
{
  int vertex = fps[i];
  P(i,0) = V(vertex,0), P(i,1) = V(vertex,1), P(i,2) = V(vertex,2);
}

// start the geodesic iso-curve signature
  int iso_start = fps[0];
  int k_iso = 20;
  double k_unit = max_geodesic/k_iso;
  std::vector<double> iso_curve_bin(k_iso, 0);
  std::vector<Eigen::VectorXd> iso_points_1, iso_points_2;
  dij_out temp_dij = dijkstra_fibo_heap(my_pair, iso_start, 0);
  for (size_t i = 0; i < m_f; i++) {
    int a = F(i,0), b = F(i,1), c = F(i,2);
    double k_a = temp_dij.dist[a];
    double k_b = temp_dij.dist[b];
    double k_c = temp_dij.dist[c];
    int radii_number = round((k_a + k_b + k_c)/(3*k_unit));
    double r_i = k_unit*radii_number;
    // std::cout << "The iso-curve bin number is " << radii_number << '\n';
    if((k_a < r_i && k_b > r_i && k_c > r_i) || (k_a > r_i && k_b < r_i && k_c < r_i))
    {
      double a1, a2;
      Eigen::VectorXd p1, p2;
      a1 = abs(r_i - k_a)/abs(k_b - k_a);
      a2 = abs(r_i - k_a)/abs(k_c - k_a);
      p1 = (1 - a1)*V.row(a) + a1*V.row(b);
      p2 = (1 - a2)*V.row(a) + a2*V.row(c);
      iso_curve_bin[radii_number] += (p1 - p2).norm();
      iso_points_1.push_back(p1), iso_points_2.push_back(p2);
    }
    else if((k_b < r_i && k_a > r_i && k_c > r_i) || (k_b > r_i && k_a < r_i && k_c < r_i))
    {
      double a1, a2;
      Eigen::VectorXd p1, p2;
      a1 = abs(r_i - k_b)/abs(k_a - k_b);
      a2 = abs(r_i - k_b)/abs(k_c - k_b);
      p1 = (1 - a1)*V.row(b) + a1*V.row(a);
      p2 = (1 - a2)*V.row(b) + a2*V.row(c);
      iso_curve_bin[radii_number] += (p1 - p2).norm();
      iso_points_1.push_back(p1), iso_points_2.push_back(p2);
    }
    else if((k_c < r_i && k_a > r_i && k_b > r_i) || (k_c > r_i && k_a < r_i && k_b < r_i))
    {
      double a1, a2;
      Eigen::VectorXd p1, p2;
      a1 = abs(r_i - k_c)/abs(k_a - k_c);
      a2 = abs(r_i - k_c)/abs(k_b - k_c);
      p1 = (1 - a1)*V.row(c) + a1*V.row(a);
      p2 = (1 - a2)*V.row(c) + a2*V.row(b);
      iso_curve_bin[radii_number] += (p1 - p2).norm();
      iso_points_1.push_back(p1), iso_points_2.push_back(p2);
    }
  }

  std::cout << "ISO-Curve Signature: " << '\n';
  for (size_t i = 0; i < iso_curve_bin.size(); i++) {
    std::cout << iso_curve_bin[i] << " ";
  }
  std::cout << '\n';

  std::ofstream myfile;
  myfile.open ("iso_and_bil_des.txt");

  for (size_t i = 0; i < fps.size(); i++) {
    std::cout << "ISO-Curve Signature" << i << ": " << "\n";
    std::vector<double> iso_curve_2 = iso_curve_signature(V, F, my_pair, fps, i);
    for (size_t j = 0; j < iso_curve_2.size(); j++) {
      std::cout << iso_curve_2[j] << " ";
    }
    std::cout << '\n';

    myfile << "ISO-Curve Signature " << i << ": \n";
    for (size_t i = 0; i < iso_curve_2.size(); i++) {
      myfile << iso_curve_2[i] << " ";
    }
    myfile << "\n";
  }




  Eigen::MatrixXd P1_iso(iso_points_1.size(), 3);
  Eigen::MatrixXd P2_iso(iso_points_2.size(), 3);

  if(iso_points_1.size() == iso_points_2.size())
  {
    std::cout << "p1 and p2 has the same length lines are drawing" << '\n';
    for (size_t i = 0; i < iso_points_1.size(); i++) {
      P1_iso.row(i) = iso_points_1[i];
      P2_iso.row(i) = iso_points_2[i];
    }
  }
  else
    std::cout << "p1 and p2 has different sizes, SHAME !!" << '\n';


  // Initialize the bilateral_hist with zeros
  std::vector<double> bilateral_hist(k_iso, 0);
  // std::cout << "Min path distance is " << my_out.min_dist << '\n';
  // std::vector<int> roi_elementes;
  std::vector<bool> is_roi_element(m_v, false);
  size_t count_roi = 0;
  double p_q_distance = my_out.min_dist;
  double distance_threshold = 0.8*my_out.min_dist;
  double bilateral_unit = p_q_distance/20;
  for (size_t i = 0; i < m_v; i++) {
    double min_geo_to_path = INT_MAX;
    for (size_t j = 0; j < parent_m.size(); j++) {
      int k = parent_m[j];
      if(path_geodesic(j, i) < min_geo_to_path && i!=k)
      {
        min_geo_to_path = path_geodesic(j, i);
      }
    }
    if(min_geo_to_path < distance_threshold &&
      path_geodesic(0, i) < p_q_distance && path_geodesic(parent_m.size()-1, i) < p_q_distance)
    {
      // roi_elementes.push_back(i);
      is_roi_element[i] = true;
      count_roi++;
      // int k_class = round(min_geo_to_path/bilateral_unit);
    }
  }

  std::cout << "ROI has " << count_roi << " elements" << '\n';
  // iterate over the faces and compute bilateral histogram
  Eigen::VectorXd Z = Eigen::VectorXd::Zero(m_f);
  for (size_t i = 0; i < m_f; i++) {
    int a = F(i,0), b = F(i,1), c = F(i,2);
    if (is_roi_element[a] && is_roi_element[b] && is_roi_element[c]) {
      double middle_bi = (path_geodesic(0, a) + path_geodesic(0, b)
        + path_geodesic(0, c))/3;
      int k_class = round(middle_bi/bilateral_unit);
      double a_b = l2_norm(V, a, b);
      double b_c = l2_norm(V, b, c);
      double a_c = l2_norm(V, a, c);
      double face_area = triangle_area_calculater(a_b, a_b, a_b);
      bilateral_hist[k_class] += face_area;
      Z(i) = 1 - middle_bi/p_q_distance;
      // std::cout << "Z value is " << Z(i) << '\n';
    }
  }



  for (size_t i = 0; i < fps.size(); i++) {
    std::cout << "ISO-Curve Signature" << i << ": " << '\n';
    std::vector<double> bilateral_hist_2 = bilateral_descriptor(V, F, my_pair, fps, i, fin);
    for (size_t j = 0; j < bilateral_hist_2.size(); j++) {
      std::cout << bilateral_hist_2[j] << " ";
    }
    std::cout << '\n';

    myfile << "bilateral histogram " << i << ": \n";
    for (size_t i = 0; i < bilateral_hist_2.size(); i++) {
      myfile << bilateral_hist_2[i] << " ";
    }
    myfile << "\n";
  }

  myfile.close();

  igl::jet(Z,true,C);

  std::cout << "ISO-K classes" << "\n";
  for (size_t i = 0; i < 20; i++) {
    std::cout << bilateral_hist[i] << " ";
  }
  std::cout << '\n';
  // // put the roi point into matrix P for visualization
  // Eigen::MatrixXd P(roi_elementes.size(), 3);
  // for(size_t i=0; i<roi_elementes.size(); i++)
  // {
  //   int vertex = roi_elementes[i];
  //   P(i,0) = V(vertex,0), P(i,1) = V(vertex,1), P(i,2) = V(vertex,2);
  // }


  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1,0,0));
  viewer.data().add_edges(P1_iso, P2_iso, Eigen::RowVector3d(1,1,0));
  viewer.data().add_points(P, Eigen::RowVector3d(0,0,1));
  viewer.data().set_colors(C);
  viewer.launch();
}
