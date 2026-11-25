#pragma once
namespace sp
{
    struct SimulationParameters
    {
        double vx_inlet, vy_inlet, total_simulation_time, frames_per_second, kinematic_viscosity, t0, dt;
    };
}