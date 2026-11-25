#include "geometric_preprocess/types.hpp"
#include "setup/fields.hpp"
namespace is
{
    void first_order_upwind_interpolation(t::MeshInfo &mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, f::Fields &fields);
    void second_order_upwind_interpolation(t::MeshInfo &mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, f::Fields &fields);
    void second_order_upwind_interpolation_with_slope_limiter(t::MeshInfo &mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, f::Fields &fields);
    void weights_based_linear_interpolation_for_velocity_gradient(t::MeshInfo& mesh_info, t::Edges& edges, f::Fields &fields);
    void orthogonality_correction_for_velocity_gradient_interpolation(t::MeshInfo& mesh_info, t::DomainTriangles& domain_triangles, f::Fields& fields);
}