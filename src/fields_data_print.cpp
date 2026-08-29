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
    void fields_data_print(t::MeshInfo& mesh_info, t::Nodes& nodes, t::DomainTriangles& domain_triangles, f::Fields& fields, int frame_number)
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
                                     << "    <Piece NumberOfPoints=\"" << mesh_info.nodes_number << "\" NumberOfCells=\"" << mesh_info.domain_triangles_number << "\">\n";

                paraview_fields_file << "      <Points>\n"
                                     << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh_info.nodes_number; ++i)
            {
                paraview_fields_file << nodes.positions[i][0] << " " << nodes.positions[i][1] << " 0 ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "      </Points>\n";

                paraview_fields_file << "      <Cells>\n"
                                     << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                paraview_fields_file << domain_triangles.nodes_IDs[i][0] << " "
                                     << domain_triangles.nodes_IDs[i][1] << " "
                                     << domain_triangles.nodes_IDs[i][2] << " ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh_info.domain_triangles_number; ++i) {
                paraview_fields_file << (i + 1) * 3 << " ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                paraview_fields_file << 5 << " ";
            }

                paraview_fields_file << "\n        </DataArray>\n"
                                     << "      </Cells>\n";

                paraview_fields_file << "      <CellData Vectors=\"Velocity\">\n"
                                     << "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                                     << "          ";

            for (i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                paraview_fields_file << fields.domain_triangles_centroids_velocities_before[i][0] << " "
                                     << fields.domain_triangles_centroids_velocities_before[i][1] << " 0 ";
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