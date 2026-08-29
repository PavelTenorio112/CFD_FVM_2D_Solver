#pragma once
namespace sp
{
    struct SimulationParameters
    {
        double vx_inlet, vy_inlet, CFL, total_simulation_time, frames_per_second, time_step;
    };
}