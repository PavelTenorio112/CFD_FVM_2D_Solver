#include "include/setup/simulation_parameters_reader.hpp"
#include "include/simulation_parameters.hpp"
#include "include/setup/fields.hpp"
#include <fstream>
#include <string>

namespace spr
{
    void simulation_parameters_reader(sp::SimulationParameters &simulation_parameters)
    {
        std::ifstream simulation_parameters_file("input/data.sp");
        std::string garbage;
        simulation_parameters_file >> garbage >> garbage 
                                   >> garbage >> simulation_parameters.vx_inlet 
                                   >> garbage >> simulation_parameters.vy_inlet
                                   >> garbage >> garbage
                                   >> garbage >> simulation_parameters.total_simulation_time
                                   >> garbage >> simulation_parameters.time_step
                                   >> garbage >> simulation_parameters.frames_per_second;
        return;
    }
}
