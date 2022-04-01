#include <iostream>
#include <cstdlib>
#include <stdio.h> 
#include "math.h"
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #Vx3
Eigen::MatrixXd V;
// Face array, #Fx3
Eigen::MatrixXi F;
Eigen::MatrixXd NF;
// Per-vertex uniform normal array, #Vx3
Eigen::MatrixXd N_uniform;
// Per-vertex area-weighted normal array, #Vx3
Eigen::MatrixXd N_area;
// Per-vertex mean-curvature normal array, #Vx3
Eigen::MatrixXd N_meanCurvature;
// Per-vertex PCA normal array, #Vx3
Eigen::MatrixXd N_PCA;
// Per-vertex quadratic fitted normal array, #Vx3
Eigen::MatrixXd N_quadraticFit;

// Per-vertex mean curvature, #Vx3
Eigen::VectorXd K_mean;
// Per-vertex Gaussian curvature, #Vx3
Eigen::VectorXd K_Gaussian;
// Per-vertex minimal principal curvature, #Vx3
Eigen::VectorXd K_min_principal;
// Per-vertex maximal principal curvature, #Vx3
Eigen::VectorXd K_max_principal;
// Per-vertex color array, #Vx3
Eigen::MatrixXd colors_per_vertex;

// Explicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_expLap;
// Bilateral smoothed vertex array, #Vx3
Eigen::MatrixXd V_bilateral;



void get_face_normals(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &Face_Normals){
    for(size_t f = 0; f < F.rows(); f++)
        {
            const Eigen::RowVector3d p0 = V.row(F(f,0));
            const Eigen::RowVector3d p1 = V.row(F(f,1));
            const Eigen::RowVector3d p2 = V.row(F(f,2));
            const Eigen::RowVector3d n0 = (p1 - p0).cross(p2 - p0);
            const Eigen::RowVector3d n1 = (p2 - p1).cross(p0 - p1);
            const Eigen::RowVector3d n2 = (p0 - p2).cross(p1 - p2);

            // careful sum
            for(int d = 0;d<3;d++)
            {
                // NF(f,d) = n0(d) + n1(d) + n2(d);
                Face_Normals(f,d) = n0(d);
            }
            Face_Normals.row(f) /= Face_Normals.row(f).norm();
        } 
}

void get_uniform_normals(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &N_uniform){
    get_face_normals(V, F, NF);
    N_uniform.setZero(V.rows(), 3);
        for(int i = 0; i< F.rows();i++)
        {
            // throw normal at each corner
            for(int j = 0; j < 3;j++)
            {
                N_uniform.row(F(i,j)) = N_uniform.row(F(i,j)) + NF.row(i);
            }
        }
        N_uniform.rowwise().normalize();
}

void get_area_normals(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &N_area){
    N_area.setZero(V.rows(), 3);
    get_face_normals(V, F, NF);
        for(int f = 0; f< F.rows();f++)
        {
            const Eigen::RowVector3d p0 = V.row(F(f,0));
            const Eigen::RowVector3d p1 = V.row(F(f,1));
            const Eigen::RowVector3d p2 = V.row(F(f,2));
            double area_face = ((p1 - p0).cross(p2 - p0)).squaredNorm()*0.5;
            // throw normal at each corner
            for(int v = 0; v < 3;v++)
            {
                N_area.row(F(f,v)) = N_area.row(F(f,v)) + area_face * NF.row(f);
            }
        }
        N_area.rowwise().normalize();
}


// void get_mean_curvature_norm(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<std::vector<int>>&Adj_list, Eigen::MatrixXd &N_area){
//     for(int v = 0; v < V.rows(); v++){ //iterate over vertices
//         std::vector<int> adjacent_list = Adj_list.at(v);
//         int ind = 0;
//         for(int ind_vertex=0;  ind_vertex < adjacent_list.size(); ind_vertex++){
//             Eigen::Vector3d current_vertex = V.row(v);
//             Eigen::Vector3d second_vertex = V.row(ind_vertex);
//             Eigen::Vector3d third_vertex = V.row((ind_vertex + 1) % adjacent_list.size());
            
//         }
//     }

// }


void get_adjacent_list(std::vector<std::vector<int>>& A, Eigen::MatrixXi &F){
      // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      // Get indices of edge: s --> d
      int s = F(i,j);
      int d = F(i,(j+1)%F.cols());
      A.at(s).push_back(d);
      A.at(d).push_back(s);
    }
  }
  
  // Remove duplicates
  for(int i=0; i<(int)A.size();++i)
  {
    std::sort(A[i].begin(), A[i].end());
    A[i].erase(std::unique(A[i].begin(), A[i].end()), A[i].end());
  }

  std::vector<std::vector<std::vector<int> > > SR; 
    SR.resize(A.size());
    
    for(int i = 0;i<F.rows();i++)
    {
      // Loop over this face
      for(int j = 0;j<F.cols();j++)
      {
        // Get indices of edge: s --> d
        int s = F(i,j);
        int d = F(i,(j+1)%F.cols());
        // Get index of opposing vertex v
        int v = F(i,(j+2)%F.cols());
        
        std::vector<int> e(2);
        e[0] = d;
        e[1] = v;
        SR[s].push_back(e);
      }
    }
    
    for(int v=0; v<(int)SR.size();++v)
    {
      std::vector<int>& vv = A.at(v);
      std::vector<std::vector<int> >& sr = SR[v];
      
      std::vector<std::vector<int> > pn = sr;
      
      // Compute previous/next for every element in sr
      for(int i=0;i<(int)sr.size();++i)
      {
        int a = sr[i][0];
        int b = sr[i][1];
        
        // search for previous
        int p = -1;
        for(int j=0;j<(int)sr.size();++j)
          if(sr[j][1] == a)
            p = j;
        pn[i][0] = p;
        
        // search for next
        int n = -1;
        for(int j=0;j<(int)sr.size();++j)
          if(sr[j][0] == b)
            n = j;
        pn[i][1] = n;
        
      }
      
      // assume manifoldness (look for beginning of a single chain)
      int c = 0;
      for(int j=0; j<=(int)sr.size();++j)
        if (pn[c][0] != -1)
          c = pn[c][0];
      
      if (pn[c][0] == -1) // border case
      {
        // finally produce the new vv relation
        for(int j=0; j<(int)sr.size();++j)
        {
          vv[j] = sr[c][0];
          if (pn[c][1] != -1)
            c = pn[c][1];
        }
        vv.back() = sr[c][1];
      }
      else
      {
        // finally produce the new vv relation
        for(int j=0; j<(int)sr.size();++j)
        {
          vv[j] = sr[c][0];
          
          c = pn[c][1];
        }
      }
    }
  }



std::vector<std::vector<int>> map_vertex_incident_faces(Eigen::MatrixXd &V, Eigen::MatrixXi &F){
  std::vector<std::vector<int>> map_vertex_faces; //map of Face indexes corresponding to the vertex index
  map_vertex_faces.resize(V.rows());
  for (int i = 0; i < F.rows(); i++){
    for (int j = 0; j < 3; j++){
      map_vertex_faces.at(F(i, j)).push_back(i);
    }
  }
  return map_vertex_faces;
}

// double area_face(int ind_face, Eigen::MatrixXi &F){



// }

double A_ij(int i, Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<std::vector<int>>&map_vertex_incident_faces){
  //i - index of the vertex
  double area_ij = 0;
  for (auto ind_face : map_vertex_incident_faces.at(i)){
    const Eigen::RowVector3d p0 = V.row(F(ind_face,0));
    const Eigen::RowVector3d p1 = V.row(F(ind_face,1));
    const Eigen::RowVector3d p2 = V.row(F(ind_face,2));
    area_ij += 1.0 / 3.0 * ((p1 - p0).cross(p2 - p0)).squaredNorm()*0.5;
  }
  return area_ij; 
}

double cot_ij(int i, int j,int ind_face, Eigen::MatrixXd &V, Eigen::MatrixXi &F){

}

Eigen::MatrixXd build_squared_distances(Eigen::MatrixXd &V, Eigen::MatrixXi &F){
  Eigen::MatrixXd squared_distances;
  squared_distances.resize(F.rows(), 3);
  for (int i = 0; i < F.rows(); i++){
    for(int j = 0; j < 3; ++j){
      squared_distances(i,0) = (V.row(F(i,1))-V.row(F(i,2))).squaredNorm(); //distance between 1 and 2
      squared_distances(i,1) = (V.row(F(i,2))-V.row(F(i,0))).squaredNorm(); //distance between 0 and 2
      squared_distances(i,2) = (V.row(F(i,0))-V.row(F(i,1))).squaredNorm(); //distance between 0 and 1
    }
  }
  return squared_distances;
}

Eigen::MatrixXd build_distances(Eigen::MatrixXd &V, Eigen::MatrixXi &F){
  Eigen::MatrixXd distances;
  distances.resize(F.rows(), 3);
  for (int i = 0; i < F.rows(); i++){
    for(int j = 0; j < 3; ++j){
      distances(i,0) = (V.row(F(i,1))-V.row(F(i,2))).norm(); //distance between 1 and 2
      distances(i,1) = (V.row(F(i,2))-V.row(F(i,0))).norm(); //distance between 0 and 2
      distances(i,2) = (V.row(F(i,0))-V.row(F(i,1))).norm(); //distance between 0 and 1
    }
  }
  return distances;
}

double cot(double a, double b, double c){
  //a and b are two incident  edges of a triangle to the angle of cot we calculate
  double cos = (a*a + b*b - c*c)/(2*a*b);
  double sin = std::sqrt(1 - cos*cos);
  double cot = cos / sin;
  return cot;
}

void get_mean_curvature_norm(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &N_meanCurvature){
  N_meanCurvature.setZero(V.rows(), 3);
  //declare size of matrix
  Eigen::SparseMatrix<double> L(V.rows(), V.rows()); //L contains all neigbouring vertex indexes j for every vertex i with weights w_{ij}
  // Eigen::SparseMatrix<double> L(V.rows(), V.rows());//declare list of non-zero elements (row, column, value)std::vector<Eigen::Triplet<double> > tripletList;   for(int= 0; i< n_faces; ++i){for(intj = 0; j < 3; ++j){doublewij= cotanWeight(i, j); //laplacianweight 1 / (2 A_i) \cot(\alpha_j)int j1= (j+1) % 3;int j2= (j+2) % 3;tripletList.push_back(Eigen::Triplet<double>(F(i, j1)), F(i, j2, wij)); tripletList.push_back(Eigen::Triplet<double>(F(i, j2)), F(i, j1, wij)); tripletList.push_back(Eigen::Triplet<double>(F(i, j1)), F(i, j1, -wij)); tripletList.push_back(Eigen::Triplet<double>(F(i, j2)), F(i, j2, -wij));  }}//construct matrix from the listL.setFromTriplets(tripletList.begin(), tripletList.end());
  //declare list of non-zero elements (row, column, value)
  std::vector<std::vector<int>> map_vertex_faces = map_vertex_incident_faces(V, F);
  // Eigen::MatrixXd squared_distances = build_squared_distances(V, F);
  Eigen::MatrixXd distances = build_distances(V, F); //build the matrix of euclidean distances
  std::vector<Eigen::Triplet<double> > tripletList;  //triplet list for all [vertex, vertex_neighbout, w_ij]
  for (int i = 0; i < F.rows(); i++){ //iterate over faces
    for(int j = 0; j < 3; ++j){ //iterate over vertexes in the face (triangle)
      double A_i = A_ij(F(i, j), V, F, map_vertex_faces);
      int j1= (j+1) % 3; //index of adjacent vertex 1
      int j2= (j+2) % 3; //index of adjacent vertex 2

      double cot_i_j1 = cot(distances(i, 0), distances(i, 2), distances(i, 1));
      double cot_i_j2 = cot(distances(i, 0), distances(i, 1), distances(i, 2));
      tripletList.push_back(Eigen::Triplet<double>(F(i, j), F(i, j1), 0.5 * cot_i_j1 / A_i)); //put the data for index F(i, j) - its 2 adjacent edges in the face i and corresponding cot_{ij}
      tripletList.push_back(Eigen::Triplet<double>(F(i, j), F(i, j2), 0.5 * cot_i_j2 / A_i));
    }
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());

  for (int k=0; k<L.outerSize(); ++k){
  for (Eigen::SparseMatrix<double>::InnerIterator it(L,k); it; ++it)
  {
    N_meanCurvature.row(it.row()) += it.value() * (V.row(it.row()) - V.row(it.col()));
  }
}
}


bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing uniform vertex normals here:
        // store in N_uniform

        NF.resize(F.rows(),3);
        get_uniform_normals(V, F, N_uniform);
        // igl::per_vertex_normals(V, F, N_uniform);
        // Set the viewer normals.
        viewer.data().set_normals(N_uniform);
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing area-weighted vertex normals here:
        // store in N_area

        // Set the viewer normals.

        // igl::per_vertex_normals(V, F, N_area);
        get_area_normals(V, F, N_area);
        viewer.data().set_normals(N_area);
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing mean-curvature vertex normals here:
        // store in N_meanCurvature
        // std::vector<std::vector<int>> A;
        // A.clear();
        // A.resize(F.maxCoeff() + 1);
        // get_adjacent_list(A, F);
        // int i = 0;
        // for (auto index_v : A){
        //     std::cout << i << std::endl;
        //     for(auto index_adj : index_v){
        //         std::cout << index_adj << ' ' ;
        //     }
        //     i += 1;
        //     std::cout<<std::endl;
        // }

        // Set the viewer normals.
        get_mean_curvature_norm(V, F, N_meanCurvature);
        viewer.data().set_normals(N_meanCurvature);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing PCA vertex normals here:
        // store in N_PCA

        // Set the viewer normals.
        viewer.data().set_normals(N_PCA);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing quadratic fitted vertex normals here:
        // store in N_quadraticFit

        // Set the viewer normals.
        viewer.data().set_normals(N_quadraticFit);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete mean curvature
        // store in K_mean

        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '7') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete Gaussian curvature
        // store in K_Gaussian

        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '8') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete minimal principal curvature
        // store in K_min_principal

        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '8') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete maximal principal curvature
        // store in K_min_principal

        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == 'E') {
        // Add your code for computing explicit Laplacian smoothing here:
        // store the smoothed vertices in V_expLap

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_expLap, F);
    }

    if (key == 'B') {
        // Add your code for computing bilateral smoothing here:
        // store the smoothed vertices in V_bilateral

        viewer.data().clear();
        viewer.data().set_mesh(V_bilateral, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOBJ(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.data().compute_normals();
  viewer.core.align_camera_center(V, F);
  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]);
    }
    else {
        filename = std::string("../data/cow.obj");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    viewer.launch();
}
