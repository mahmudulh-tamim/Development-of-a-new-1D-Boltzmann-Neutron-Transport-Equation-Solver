function flux_new=group_flux(flux_old,k_old)

%spatial discretization

slab_length=2;
mesh_number=100;
del_x=2/mesh_number;
x=(0:del_x:2)';
edge_count=length(x);
mesh_count=edge_count-1; % it should be equal to mesh_number

mesh_length=zeros(mesh_count,1);
mesh_length(1:end,1)=del_x;

flux_new=zeros(mesh_count,3);

for g=1:3
    flux_new(:,g)=source_iteration(flux_old,k_old,g);
end
