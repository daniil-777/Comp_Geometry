#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <cmath>
#include <igl/slice.h>

#define MAX 100000

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;


class Shape_Handler
{
public:
    Shape_Handler();
    void init_global_variable();
    void createGrid();
    void evaluateImplicitFunc();
    void getLines();
    // bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
    void neighbors(Eigen::RowVector3d p);
    int coordinates_to_spatial_idx(int x, int y, int z);
    int grid_coordinates_to_idx(int x, int y, int z);
    Eigen::VectorXd wendLand(std::vector<double> r);
    int Grid_query_index(Eigen::RowVector3d p);
    void createGrid_query();
    std::vector<int> point_to_grid(Eigen::RowVector3d p);
    int idx_closest_point(Eigen::RowVector3d p);
    bool callback_load_mesh(Viewer& viewer,string filename);
    // void neighbors(Eigen::RowVector3d p);
    Eigen::VectorXd get_polynome_coef(Eigen::MatrixXd data_array, int index, int degree);
    


// protected:

// 	Eigen::MatrixXd wireframe_edges_from;
// 	Eigen::MatrixXd wireframe_edges_to;
// 	vector<Eigen::Vector3d> m_grid_matrix_positions;
//     Eigen::RowVector3d m_wireframe_color;


};