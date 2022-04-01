#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
<<<<<<< HEAD
#include <cmath>
#include <igl/slice.h>
#include "main.h"
#define MAX 100000
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
<<<<<<< HEAD
Eigen::MatrixXd P; //with general index
=======
Eigen::MatrixXd P;
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;
<<<<<<< HEAD
double diag;
double stepRate = 0.1;
double scale = 1.2;
double grid_step_x, grid_step_y, grid_step_z;
double offset_view = 1.1;
double epsConst = 0.05;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;
double wendlandRadiusDiag;

// Parameter: grid resolution
int resolutionX = 20;
int resolutionY = 20;
int resolutionZ = 20;



// Intermediate result: grid points, at which the implicit function will be evaluated, #G x3
=======

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
int resolution = 20;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

<<<<<<< HEAD
//new structure for grid
std::vector<std::vector<int>> Grid_query;
int dim_Grid_query_x, dim_Grid_query_y, dim_Grid_query_z;
double Grid_query_step;
bool Normal_Constraints = false;
bool PCA = false; //using PCA
bool BRUTE_FORCE = false; //using simple iteration over all points
double Min_Distance_Init= 1000000.0;
Eigen::VectorXd saveConstrValues;
Eigen::MatrixXd closepoints;
Eigen::VectorXi neighbors_points;
string shapeName;
Eigen::MatrixXd tempP, tempN;
Eigen::Matrix3d temp;
// Global variable for grid bounds
Eigen::RowVector3d  bb_min, bb_max, dim;
bool singleton_PCA = true;
std::vector<double> distanceVector;



void init_global_variable(){

    tempP = P;
    tempN = N;

    if(PCA && singleton_PCA) {
        Eigen::MatrixXd centerPoints = P.rowwise() - P.colwise().mean();
        Eigen::MatrixXd cov = centerPoints.adjoint() * centerPoints;
        //cov = cov / (P.rows() - 1);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(cov);
        Eigen::MatrixXd eigenVectors = eigen.eigenvectors();
        eigenVectors = eigenVectors.rightCols(3);
        tempP = (P * eigenVectors);
        tempN = (N * eigenVectors);
        singleton_PCA = false;
    }

    P = tempP;
    N = tempN;

    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();
    dim = (bb_max-bb_min)*scale;
    diag = dim.norm();
    wendlandRadiusDiag = wendlandRadius*diag;
}
int coordinates_to_spatial_idx(int x, int y, int z)
{
    return x + dim_Grid_query_x * y + dim_Grid_query_x * dim_Grid_query_y * z ;
}
int grid_coordinates_to_idx(int x, int y, int z)
{
    return x + resolutionX*y + resolutionX*resolutionY*z;
}

//returns the wendLand function
Eigen::VectorXd wendLand(std::vector<double> r){
    Eigen::VectorXd dist = Eigen::VectorXd::Map(&r[0], r.size());
    Eigen::VectorXd result = (1 - (dist.array()/(wendlandRadiusDiag))).pow(4) * (4 * dist.array()/(wendlandRadiusDiag) + 1);
    return result;
}

void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);


    //Global variable
    grid_step_x = dim[0] / (double)(resolutionX - 1);
    grid_step_y = dim[1] / (double)(resolutionY - 1);
    grid_step_z = dim[2] / (double)(resolutionZ - 1);

    //enlarge the grid
    Eigen::RowVector3d diff(3);
    diff << offset_view * grid_step_x, offset_view * grid_step_y, offset_view * grid_step_z;

    grid_points.resize(resolutionX*resolutionY*resolutionZ, 3);
    for (int i = 0; i < resolutionX; ++i) {
        for (int j = 0; j < resolutionY; ++j) {
            for (int k = 0; k < resolutionZ; ++k) {
                grid_points.row(grid_coordinates_to_idx(i, j, k)) = bb_min -diff + Eigen::RowVector3d(i * grid_step_x, j * grid_step_y, k * grid_step_z);
            }
        }
    }

}

int Grid_query_index(Eigen::RowVector3d p){
    //change ccordinated of a point to "disctrete coordinates"
    Eigen::RowVector3d p_discrete = (p - bb_min)/Grid_query_step; //coordinates in the system with bb_min as a start and bb_min as a step
    int X = floor(p_discrete.x());
    int Y = floor(p_discrete.y());
    int Z = floor((p_discrete.z()));
    return coordinates_to_spatial_idx(X, Y, Z);
}

void createGrid_query(){
    Grid_query.clear();
    Grid_query_step = stepRate * diag;

    dim_Grid_query_x = ceil(dim.x() / Grid_query_step);
    dim_Grid_query_y = ceil(dim.y() / Grid_query_step);
    dim_Grid_query_z = ceil(dim.z() / Grid_query_step);
    int Grid_query_size = dim_Grid_query_x * dim_Grid_query_y * dim_Grid_query_z;
    Grid_query.resize(Grid_query_size);

    for (int i = 0; i < P.rows(); ++i) {
        int spatial_index = Grid_query_index(P.row(i));
        // move index of the point in the Points array to the slot with "spatial_index" in Grid_query
        Grid_query[spatial_index].push_back(i); 
    }
}

std::vector<int> point_to_grid(Eigen::RowVector3d p){
    Eigen::Vector3d p_discrete = (p - bb_min)/Grid_query_step;
    int X = floor(p_discrete[0]);
    int Y = floor(p_discrete[1]);
    int Z = floor(p_discrete[2]);
    std::vector<int> result{X,Y,Z};
    return result;
}

int idx_closest_point(Eigen::RowVector3d p){
    //find the general index of the point closest to p
    double minimumDistance = Min_Distance_Init;
    std::vector<int> xyz = point_to_grid(p);
    int X = xyz[0];
    int Y = xyz[1];
    int Z = xyz[2];
    int distance_search = 1;
    //find coordinates of neighbouring grid cells around point p
    int minQuery_X = max(0, X-distance_search);  //"left" X
    int maxQuery_X = min(X + 1 + distance_search, dim_Grid_query_x); //"right" X
    int minQuery_Y = max(0, Y-distance_search);
    int maxQuery_Y = min(Y + 1 + distance_search, dim_Grid_query_y);
    int minQuery_Z = max(0, Z-distance_search);
    int maxQuery_Z = min(Z + 1 + distance_search, dim_Grid_query_z);
    int minDistanceIndex = -1;

    if(!BRUTE_FORCE){
    for (int i = minQuery_X; i < maxQuery_X ; ++i) {
        for (int j = minQuery_Y; j < maxQuery_Y; ++j) {
            for (int k = minQuery_Z; k < maxQuery_Z ; ++k) {
                int spatial_index = coordinates_to_spatial_idx(i, j, k); // get global (flattened index for the spatial grid)
                //iterate over points 
                for (int u = 0; u < Grid_query[spatial_index].size(); ++u) {
                    Eigen::RowVector3d dist = (p - P.row(Grid_query[spatial_index][u]));
                    if (dist.norm() < minimumDistance){
                        minDistanceIndex = Grid_query[spatial_index][u];
                        minimumDistance = dist.norm();
                    }
                }
            }
        }
    }
    minimumDistance = Min_Distance_Init;
    return minDistanceIndex;
}else {
    for (int i = 0; i <  P.rows(); ++i) {
        Eigen::RowVector3d dist = (p - P.row(i));
        if (dist.norm() < minimumDistance){
            minDistanceIndex = i;
            minimumDistance = dist.norm();
        }
    }
    }
    minimumDistance = Min_Distance_Init;
    return minDistanceIndex;
}


void neighbors(Eigen::RowVector3d p){
    if(!Normal_Constraints){
        saveConstrValues.setZero(constrained_points.rows(), 1);
        neighbors_points.setZero(constrained_points.rows());
        distanceVector.clear();
        int cont = 0;
        if(!BRUTE_FORCE){
        std::vector<int> xyz = point_to_grid(p);
        int X = xyz[0];
        int Y = xyz[1];
        int Z = xyz[2];
        int distance_search = floor(wendlandRadiusDiag/Grid_query_step)+1;
        int minQuery_X = max(0, X-distance_search);
        int maxQuery_X = min(X + 1 + distance_search, dim_Grid_query_x);
        int minQuery_Y = max(0, Y-distance_search);
        int maxQuery_Y = min(Y + 1 + distance_search, dim_Grid_query_y);
        int minQuery_Z = max(0, Z-distance_search);
        int maxQuery_Z = min(Z + 1 + distance_search, dim_Grid_query_z);
        
        
        
        for (int i = minQuery_X; i < maxQuery_X ; ++i) {
            for (int j = minQuery_Y; j < maxQuery_Y; ++j) {
                for (int k = minQuery_Z; k < maxQuery_Z ; ++k) {
                    int spatial_index = coordinates_to_spatial_idx(i, j, k);
                    for (int u = 0; u < Grid_query[spatial_index].size(); ++u) {
                        double dist = (p - constrained_points.row(Grid_query[spatial_index][u])).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = Grid_query[spatial_index][u];
                            saveConstrValues.row(cont) = constrained_values.row(Grid_query[spatial_index][u]);
                            cont++;
                        }
                        dist = (p - constrained_points.row(P.rows() + Grid_query[spatial_index][u])).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = P.rows() + Grid_query[spatial_index][u];
                            saveConstrValues.row(cont) = constrained_values.row(P.rows() + Grid_query[spatial_index][u]);
                            cont++;
                        }
                        dist = (p - constrained_points.row(2*P.rows() + Grid_query[spatial_index][u])).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = 2*P.rows() + Grid_query[spatial_index][u];
                            saveConstrValues.row(cont) = constrained_values.row(2*P.rows() + Grid_query[spatial_index][u]);
                            cont++;
                        }
                    }
                }
            }
        }
        neighbors_points.conservativeResize(cont);
        saveConstrValues.conservativeResize(cont, 1);
        } else {
            for (int i = 0; i < P.rows(); ++i) {
                        double dist = (p - constrained_points.row(i)).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = i;
                            saveConstrValues.row(cont) = constrained_values.row(i);
                            cont++;
                        }
                        dist = (p - constrained_points.row(P.rows() + i)).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = P.rows() + i;
                            saveConstrValues.row(cont) = constrained_values.row(P.rows() + i);
                            cont++;
                        }
                        dist = (p - constrained_points.row(2*P.rows() + i)).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = 2*P.rows() + i;
                            saveConstrValues.row(cont) = constrained_values.row(2*P.rows() + i);
                            cont++;
                        }
   
        }
        neighbors_points.conservativeResize(cont);
        saveConstrValues.conservativeResize(cont, 1);
        }
    }
    else{
        saveConstrValues.setZero(P.rows(), 1);
        neighbors_points.setZero(P.rows());
        distanceVector.clear();
        
        std::vector<int> xyz = point_to_grid(p);
        int X = xyz[0];
        int Y = xyz[1];
        int Z = xyz[2];
        int distance_search = floor(wendlandRadiusDiag/Grid_query_step)+1;
        int minQuery_X = max(0, X-distance_search);
        int maxQuery_X = min(X + 1 + distance_search, dim_Grid_query_x);
        int minQuery_Y = max(0, Y-distance_search);
        int maxQuery_Y = min(Y + 1 + distance_search, dim_Grid_query_y);
        int minQuery_Z = max(0, Z-distance_search);
        int maxQuery_Z = min(Z + 1 + distance_search, dim_Grid_query_z);

        int cont = 0;
        for (int i = minQuery_X; i < maxQuery_X ; ++i) {
            for (int j = minQuery_Y; j < maxQuery_Y; ++j) {
                for (int k = minQuery_Z; k < maxQuery_Z ; ++k) {
                    int spatial_index = coordinates_to_spatial_idx(i, j, k);
                    for (int u = 0; u < Grid_query[spatial_index].size(); ++u) {
                        double dist = (p - P.row(Grid_query[spatial_index][u])).norm();
                        if (dist < wendlandRadiusDiag){
                            distanceVector.push_back(dist);
                            neighbors_points(cont) = Grid_query[spatial_index][u];
                            saveConstrValues.row(cont) = constrained_values.row(Grid_query[spatial_index][u]); 
                            cont++;
                        }
                    }
                }
            }
        }
        neighbors_points.conservativeResize(cont);
        saveConstrValues.conservativeResize(cont, 1);
    }

}

Eigen::VectorXd get_polynome_coef(Eigen::MatrixXd data_array, int index, int degree){
    //get basis functions for the array of indicies for the point cloud, index and degree of polymomial 

    Eigen::VectorXd answer(1);
    if(degree == 0){
        answer << 1;
        return answer;
    }else if(degree == 1){
        answer.resize(4);
        answer << 1, data_array(index, 0), data_array(index, 1), data_array(index, 2);
        return answer;
    }else if(degree == 2){
        answer.resize(10);
        answer << 1, data_array(index, 0), data_array(
                                        index, 1), data_array(index, 2),
                                        pow(data_array(index, 0), 2), pow(
                                        data_array(index, 1), 2), pow(
                                        data_array(index, 2), 2),
                                        data_array(index, 0) *
                                        data_array(index, 1),
                                        data_array(index, 1) *
                                        data_array(index, 2),
                                        data_array(index, 0) *
                                        data_array(index, 2);
        return answer;
    }
}


void evaluateImplicitFunc(){

    grid_values.resize(resolutionX*resolutionY*resolutionZ);
    grid_values.setZero(resolutionX*resolutionY*resolutionZ);

    for (int i = 0; i < resolutionX; ++i) {
        for (int j = 0; j < resolutionY; ++j) {
            for (int k = 0; k < resolutionZ; ++k) {
                int grid_index = grid_coordinates_to_idx(i,j,k);

                Eigen::MatrixXd b;
                Eigen::VectorXd weightVec;


                neighbors(grid_points.row(grid_index));

                int number_neighbours = neighbors_points.size();
                if(number_neighbours == 0)
                    grid_values[grid_index] = MAX;
                else{
                    Eigen::VectorXd b_x(1);
                    weightVec = wendLand(distanceVector);
                    if(polyDegree == 0){
                        b.resize(number_neighbours, 1);
                        for (int i = 0; i < number_neighbours; i++) {
                            // b.row(i) << 1;
                            b.row(i) << get_polynome_coef(constrained_points, i, 0);
                        }
                        b_x << 1;
                    }
                    else if(polyDegree == 1){
                        b.resize(number_neighbours, 4);
                        for (int i = 0; i < number_neighbours; i++) {
                            if(!Normal_Constraints) {
                                
                                b.row(i) = get_polynome_coef(constrained_points, neighbors_points(i), 1);
                                // b.row(i) << 1, constrained_points(neighbors_points(i), 0), constrained_points(
                                //         neighbors_points(i), 1), constrained_points(neighbors_points(i), 2);
                            }
                            else{
   
                                b.row(i) << get_polynome_coef(P, neighbors_points(i), 1); //just se initial points
                            }
                        }
                        b_x.resize(4);
                        b_x << 1, grid_points(grid_index,0), grid_points(grid_index,1), grid_points(grid_index,2);
                    }
                    else if(polyDegree == 2){
                        b.resize(number_neighbours, 10);
                        for (int i = 0; i < number_neighbours; i++) {
                            if(!Normal_Constraints) {
                                
                                b.row(i) = get_polynome_coef(constrained_points, neighbors_points(i), 2);
                            }
                            else{
                                b.row(i) = get_polynome_coef(P, neighbors_points(i), 2);
                            }
                        }
                        b_x.resize(10);
                        b_x << 1, grid_points(grid_index,0), grid_points(grid_index,1), grid_points(grid_index,2),
                                pow(grid_points(grid_index,0),2),  pow(grid_points(grid_index,1),2),  pow(grid_points(grid_index,2), 2),
                                grid_points(grid_index,0)* grid_points(grid_index,1),  grid_points(grid_index,1)* grid_points(grid_index,2),  grid_points(grid_index,0)* grid_points(grid_index,2);
                    }

                    if(Normal_Constraints) {
                        for(int i=0; i<saveConstrValues.size(); i++) {
                            saveConstrValues(i) += (grid_points.row(grid_index) - P.row(neighbors_points(i))).dot(N.row(neighbors_points(i)%N.rows())); // (x - p)n
                        }
                    }

                    Eigen::MatrixXd A = weightVec.asDiagonal()*b; // W*B
                    Eigen::VectorXd c = A.colPivHouseholderQr().solve((weightVec).asDiagonal()*saveConstrValues); // (WB)^(-1)*W*phi
                    grid_values[grid_index] = b_x.dot(c);
                }
=======
// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines. resize(0, 6);
    grid_values.resize(0);
    V. resize(0, 3);
    F. resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc() {
    // Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    double radius = 0.5 * (bb_max - bb_min).minCoeff();

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                // Value at (x,y,z) = implicit function for the sphere
                grid_values[index] = (grid_points.row(index) - center).norm() - radius;
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
            }
        }
    }
}

<<<<<<< HEAD

    

=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

<<<<<<< HEAD
    for (unsigned int x = 0; x<resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionX - 1) {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1) {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1) {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionX * resolutionY;
=======
    for (unsigned int x = 0; x<resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1) {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1) {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1) {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }
<<<<<<< HEAD
=======

>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
<<<<<<< HEAD
        init_global_variable();

        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }
    if (key == '2') {
        // Show all constraints
        init_global_variable();
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        // Add your code for computing auxiliary constraint points here

        //fast data structure of query access to the point
        createGrid_query();
        // normalize normal of every point
        N.rowwise().normalize();
        int dimP = P.rows();
        //values of constrained points (inner, zero level, outer), so it is numberPoints*3
        constrained_values.resize(dimP * 3);
        //coordinates of constrained points (inner, zero level, outer), so it is numberPoints*3,3
        constrained_points.resize(dimP * 3, 3);


        double eps = epsConst * diag; //define epsilon

        //iterate over initial points P
        for (int i = 0; i < dimP; ++i) {
            constrained_points.row(i) = P.row(i);
            constrained_values(i) = 0;
            constrained_points.row(i+dimP) = P.row(i) + eps*N.row(i);
            //check whether the point in the direction of eps*N is the closest t i-th point in P point cloud
            while (idx_closest_point(P.row(i) + eps*N.row(i)) != i){
                eps *= 0.5;
            }
            //assign eps scalar value in that point
            constrained_points.row(i+dimP) = P.row(i) + eps*N.row(i);
            constrained_values(i+dimP) = eps;
            eps = epsConst * diag;
            constrained_points.row(i+2*dimP) = P.row(i) - eps*N.row(i);
            while (idx_closest_point(P.row(i) - eps*N.row(i)) != i){
                eps *= 0.5;
            }
            constrained_points.row(i+2*dimP) = P.row(i) - eps*N.row(i);
            constrained_values(i+2*dimP) = -eps;
        }
        // Display all points (zero level, outer and inner points)
        viewer.data().clear();
        viewer.core.align_camera_center(constrained_points);
        viewer.data().point_size = 6;
        //blue zero level 
        viewer.data().add_points(constrained_points.block(0, 0, dimP, 3), Eigen::RowVector3d(0,0,1));
        //red outer points
        viewer.data().add_points(constrained_points.block(dimP, 0, dimP, 3), Eigen::RowVector3d(1,0,0));
        //green inner points
        viewer.data().add_points(constrained_points.block(dimP * 2, 0, dimP, 3), Eigen::RowVector3d(0,1,0));
    }
    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
        if (constrained_points.rows() == 0)
            callback_key_down(viewer, '2', 0);
=======
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0,0,0));
    }

    if (key == '2') {
        // Show all constraints
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        // Add your code for computing auxiliary constraint points here
        // Add code for displaying all points, as above
    }

    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

<<<<<<< HEAD
        //cout << constrained_values << endl;

=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();
<<<<<<< HEAD
        //cout << "ok finito la prima stampa" << saveConstrValues << endl;
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i) {
            double value = grid_values(i);
            if (value < 0) {
                grid_colors(i, 1) = 1;
            }
            else {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }
<<<<<<< HEAD
=======

>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
<<<<<<< HEAD
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/

=======
                              grid_lines.block(0, 3, grid_lines.rows(), 3),
                              Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    }

    if (key == '4') {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
<<<<<<< HEAD
            /* cerr << "Not enough data for Marching Cubes !" << endl;
             return true;*/
            callback_key_down(viewer, '3', 0);
            viewer.data().clear();
        }


        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
=======
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        if (V.rows() == 0) {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

<<<<<<< HEAD

=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
    }

<<<<<<< HEAD
    //Save
    if (key == '0') {
        igl::writeOFF("../results/res.off", V, F);
        cout << "Saved reconstructed " << shapeName << endl;
    }
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
<<<<<<< HEAD
    igl::readOFF(filename,P,F,N);
    callback_key_down(viewer,'1',0);
    return true;
}


int main(int argc, char *argv[]){
    if (argc != 2) {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/sphere.off",P,F,N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1],P,F,N);    
    }
=======
  igl::readOFF(filename,P,F,N);
  callback_key_down(viewer,'1',0);
  return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
      cout << "Usage ex2_bin <mesh.off>" << endl;
      igl::readOFF("../data/sphere.off",P,F,N);
    }
	  else
	  {
		  // Read points and normals
		  igl::readOFF(argv[1],P,F,N);
	  }
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
<<<<<<< HEAD
    // Shape_Handler shape_handler;
    // Shape_Handler *shape_handler = NULL;
    // shape_handler = new Shape_Handler();
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
<<<<<<< HEAD
        // Draw parent menu content
        menu.draw_viewer_menu();
        // callback_key_down(viewer,'2',0);

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputInt("X resolution", &resolutionX, 0, 0);
            ImGui::InputInt("Y resolution", &resolutionY, 0, 0);
            ImGui::InputInt("Z resolution", &resolutionZ, 0, 0);
            ImGui::InputInt("Degree of Polynom", &polyDegree, 0, 0);
            ImGui::InputDouble("Wendland Radius Default", &wendlandRadius, 0, 0);
            ImGui::Checkbox("Use Brutal Force", &BRUTE_FORCE);
            ImGui::Checkbox("PCA computation", &PCA);
            ImGui::Checkbox("Normal Constraints", &Normal_Constraints);
            ImGui::InputDouble("offset view", &offset_view, 0, 0);
            ImGui::InputDouble("Epsilon", &epsConst, 0, 0);

            if (ImGui::Button("Reset Grid", ImVec2(-1,0)))
            {
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer,'3',0);
            }

            // TODO: Add more parameters to tweak here...
        }

=======
      // Draw parent menu content
      menu.draw_viewer_menu();

      // Add new group
      if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
      {
        // Expose variable directly ...
        ImGui::InputInt("Resolution", &resolution, 0, 0);
        if (ImGui::Button("Reset Grid", ImVec2(-1,0)))
        {
          std::cout << "ResetGrid\n";
          // Recreate the grid
          createGrid();
          // Switch view to show the grid
          callback_key_down(viewer,'3',0);
        }

        // TODO: Add more parameters to tweak here...
      }

>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
<<<<<<< HEAD
}
=======
}
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
