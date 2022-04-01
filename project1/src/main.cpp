#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
<<<<<<< HEAD
/*** insert any libigl headers here ***/
#include <imgui/imgui.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/flip_edge.h>

#include <iostream>
#include <fstream>
=======
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

<<<<<<< HEAD
#define PI 3.14159265
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
<<<<<<< HEAD
        ofstream myfile;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
        myfile.open("data_structure_1.txt");
        myfile << "Vertex              adjacent faces.\n";
        for (int i = 0; i < VF.size(); i++)
        {
            //I am writing inf to a text file and then use in Readme
            // cout << "Vertex " << i << " -------adjacent faces ";
            myfile << "|" << i << "|";
            // for (int j = 0; j < VF[i].size(); j++)
            for(auto j: VF[i])
            {
                cout << j << " ";
                myfile << j << "  ";
            }
            myfile << "|\n";
            
            cout << endl;
        }
        myfile.close();
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        
=======
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
<<<<<<< HEAD
        igl::adjacency_list(F, VV);
        ofstream myfile;
        myfile.open("data_structure_2.txt");
        myfile << "Vertex              adjacent vertixes.\n";
        for (int i = 0; i < VV.size(); i++)
        {
            cout << "Vertex " << i << " connected to: vertices ";
            myfile << "|" << i << "|";
            for (auto current_vertex: VV[i])
            // for (int j = 0; j < VV[i].size(); j++)

            {
                // cout << VV[i][j] << " ";
                cout << current_vertex << " ";
                myfile << current_vertex << "  ";
            }
            myfile << "|\n";
            cout << endl;
        }
        myfile.close();

=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
<<<<<<< HEAD
        igl::per_face_normals(V, F, FN);
=======

>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
<<<<<<< HEAD
   
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V, F, VN);
        // Set the viewer normals.
        viewer.data().set_normals(VN);
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
        // Add your code for computing per-vertex normals here: store in VN.

        // Set the viewer normals.
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
<<<<<<< HEAD
        const double threshold = 15;
        igl::per_corner_normals(V, F, threshold, CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
        // Add your code for computing per-corner normals here: store in CN.
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6

        //Set the viewer normals
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
<<<<<<< HEAD
        component_colors_per_face.setZero(F.rows(), 3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F, cid);
        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        igl::jet(cid, true, component_colors_per_face);
        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
        std::map<int, int> map_components_frequence;
        for (int index_face = 0; index_face < F.rows(); index_face++)
        {
            int index_component = cid[index_face];
            // std:: cout << "fi" << fi;
            map_components_frequence[index_component]++;
        }
        ofstream myfile;
        //for my Readme
        myfile.open("components_statistics.txt");
        // myfile << "%s\n";
        myfile << "number of components is: " << map_components_frequence.size() << "\n";
        myfile << "id of a component    number of faces.\n";
        
        std::cout << "number of components is: " << map_components_frequence.size() << std::endl;
        std::cout << "id of a component   number of faces" << std::endl;
        for (auto index_map = map_components_frequence.begin(); index_map != map_components_frequence.end(); ++index_map)
        // for (auto index : component_count)
        {
            // myfile.open("components_statistics.txt");
            myfile << "|" << index_map->first << " | " << index_map->second << "|" << "\n";
            std::cout << index_map->first << "\t" << "\t    \t"<< index_map->second << std::endl;
        }
        myfile.close();
    }

    if (key == '7') {
        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;
        // Add your code for sqrt(3) subdivision here.
        Eigen::MatrixXd BC;
        //creation of baricenter - step 1
        igl::barycenter(V, F, BC);

        //new set of vertexes - [V M] - step 1
        Vout.setZero(V.rows() + BC.rows(), V.cols());
        Vout << V, BC;

        Eigen::MatrixXi EV, FE, EF;
        //according to /igl/edge_topology.cpp/hpp:
        //F - list of Faces N*3
        //EV - list of pair of vertexes for every edge
        //FE - list of pairs of Faces neighboring every edge 
        igl::edge_topology(V, F, EV, FE, EF);

        //since every face is replaced with 3 faces as a result of baricenter-vertexes connection; edges flipping doesn't saves the number of faces
        Fout.setZero(F.rows() * 3, F.cols());

        igl::adjacency_list(F, VV);
        Eigen::MatrixXd sum_neighbouring_vertexes;
        sum_neighbouring_vertexes.setZero(V.rows(), V.cols());
        //2) moving + averaging
        for (int i = 0; i < V.rows(); i++)
        {
            //number of neighboring vertexes for each vertex
            int n = VV[i].size();
            // sum_neighbouring_vertexes.setZero();
            for (int j = 0; j < n; j++)
            {
                sum_neighbouring_vertexes.row(i) += V.row(VV[i][j]); //VV[i][j] is index of j-th vector adjacent to i-th vector
            }
            double a = (4 - 2 * cos((2 * M_PI) / n)) / 9;
            Vout.row(i) << (1 - a) * V.row(i) + a / n * sum_neighbouring_vertexes.row(i);
        }


        //Fout.resize(F.rows()*3,3);
        int index_new_face = 0, vertex_left, vertex_right, vertex_baricenter, baricenter_left, baricenter_right;
        for (int i = 0; i < EV.rows(); i++)
        {
            if (EF(i, 0) == -1 || EF(i, 1) == -1)
            { 
                //at the border there is only one face so we have -1
                //connect both vertices and barycenter
                vertex_left = EV(i, 0);
                vertex_right = EV(i, 1);

                //for sure there would be -1, so it would be canceled out with 1
                //I store firstly vertices (with numeration as rows of faces) and then baricenters, they are "shifted" to the number of faces below
                //to get the index of baricenter I add the index of corresponding row of faces and the whole number of faces

                vertex_baricenter = 1 + EF(i, 0) + EF(i, 1) + V.rows();
                Fout.row(index_new_face) = Eigen::Vector3i(vertex_baricenter, vertex_right, vertex_left);

                index_new_face++;
            }
            else
            { //without border
                
                //I store firstly vertices (with numeration as rows of faces) and then baricenters, they are "shifted" to the number of faces below
                //to get the index of baricenter I add the index of corresponding row of faces and the whole number of faces
                baricenter_left = EF(i, 0) + V.rows();
                baricenter_right = EF(i, 1) + V.rows();
                vertex_left = EV(i, 0);
                vertex_right = EV(i, 1);
    
                //add to a list of faces face made of both edge's vertexes and a baricenter in one side ("left")
                //we are flipping original edge so new face has one "moved vertex P and two baricenters-vertexes"
                Fout.row(index_new_face) = Eigen::Vector3i(baricenter_right, baricenter_left, vertex_left);
        
                //add to a list of faces face made of both edge's vertexes and a baricenter in another side ("right")
                Fout.row(index_new_face + 1) = Eigen::Vector3i(baricenter_left, baricenter_right, vertex_right);
                //cout << Fout.row(c + 1) << endl;
                index_new_face += 2;
            }
        }
        
        // Set up the viewer to display the new mesh
        V = Vout;
        F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);

=======
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.

        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
        Eigen::MatrixXd Vout=V;
        Eigen::MatrixXi Fout=F;
        // Add your code for sqrt(3) subdivision here.

        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    }

    return true;
}

<<<<<<< HEAD
    
=======
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
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
<<<<<<< HEAD
    // viewer.load_mesh = load_mesh;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

=======
    
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]);
    }
    else {
<<<<<<< HEAD
        filename = std::string("../data/coffecup.off");
    }
    load_mesh(viewer,filename,V,F);
    

    callback_key_down(viewer, '2', 0);
    // callback_key_down(viewer, '7', 0);
    // callback_key_down(viewer, '7', 0);
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();
    };


=======
        filename = std::string("../data/bunny.off");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
>>>>>>> 66a8c78257da8ab818e78717e519cb108ab781c6
    viewer.launch();
}
