function flux_new=group_flux(flux_old,k_old)
%spatial discretiization
%uranium

half_slab_length_U=0.341011;
mesh_number_U=100;

mesh_length_U=half_slab_length_U/mesh_number_U;
half_x_U=(0:mesh_length_U:half_slab_length_U)';

half_edge_count_U=length(half_x_U);
half_mesh_count_U=half_edge_count_U-1;

%H2O
half_slab_length_H2O=0.751023;
mesh_number_H2O=100;

mesh_length_H2O=half_slab_length_H2O/mesh_number_H2O;
half_x_H2O=(half_slab_length_U:mesh_length_H2O:half_slab_length_H2O+half_slab_length_U)';

half_edge_count_H2O=length(half_x_H2O);
half_mesh_count_H2O=half_edge_count_H2O-1;

%negative region
neg_half_x_U=-flip(half_x_U,1);
neg_half_x_H2O=-flip(half_x_H2O,1);


%spatial ordinates

x=cat(1,neg_half_x_H2O(1:end-1,1),neg_half_x_U(1:end-1,1),half_x_U(1:end-1,1),half_x_H2O(1:end,1));
edge_count=length(x);
mesh_count=edge_count-1;

mesh_length=zeros(mesh_count,1);

mesh_length(1:half_mesh_count_H2O,1)=mesh_length_H2O;
mesh_length(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,1)=mesh_length_U;
mesh_length(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,1)=mesh_length_U;
mesh_length(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,1)=mesh_length_H2O;


flux_new=zeros(mesh_count,2);

for g=1:2
    flux_new(:,g)=source_iteration(flux_old,k_old,g);
end
