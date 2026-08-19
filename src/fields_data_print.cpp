#include "include/fields_data_print.hpp"
#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace fdp
{
    void fields_data_print(t::MeshData& mesh, fds::Fields& fields, int frame_number)
    {
        int i;
        std::string file_number;
        std::ostringstream file_number_stream;
        file_number_stream << std::setw(4) << std::setfill('0') << frame_number;
        std::string file_number_string = file_number_stream.str();
        std::string file_name = "output/paraview/H" + file_number_string + ".vtu";
        std::ofstream paraview_fields_file(file_name);
    
        if(paraview_fields_file.is_open())
        {
                paraview_fields_file << "<?xml version=\"1.0\"?>\n"
                                     << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                                     << "  <UnstructuredGrid>\n"
                                     << "    <Piece NumberOfPoints=\"" << mesh.nodes_number << "\" NumberOfCells=\"" << mesh.domain_triangles_number << "\">\n";

                paraview_fields_file << "      <Points>\n"
                                     << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh.nodes_number; ++i)
            {
                paraview_fields_file << mesh.nodes[i].position[0] << " " << mesh.nodes[i].position[1] << " 0 ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "      </Points>\n";

                paraview_fields_file << "      <Cells>\n"
                                     << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh.domain_triangles_number; ++i)
            {
                paraview_fields_file << mesh.domain_triangles[i].nodes_ids[0] << " "
                                     << mesh.domain_triangles[i].nodes_ids[1] << " "
                                     << mesh.domain_triangles[i].nodes_ids[2] << " ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh.domain_triangles_number; ++i) {
                paraview_fields_file << (i + 1) * 3 << " ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh.domain_triangles_number; ++i)
            {
                paraview_fields_file << 5 << " ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "      </Cells>\n";

                paraview_fields_file << "      <CellData Vectors=\"Velocity\">\n"
                                     << "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh.domain_triangles_number; ++i) 
            {
                paraview_fields_file << fields.triangles_centroids_vx_before[i] << " "
                                     << fields.triangles_centroids_vy_before[i] << " 0 ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "      </CellData>\n";

                paraview_fields_file << "    </Piece>\n"
                                     << "  </UnstructuredGrid>\n"
                                     << "</VTKFile>\n";

            paraview_fields_file.close();
        }
        else
        {
            std::cout << "There was an error printing fields data. " << std::endl;
        }
        return;
    }
}