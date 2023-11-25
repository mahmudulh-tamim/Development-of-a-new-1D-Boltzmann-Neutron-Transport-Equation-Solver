function flux_new=source_iteration(flux_old,k_old)


tol=10^-6;

sigma_t_1=1;
sigma_a_1=0.3;
sigma_s_1=1-sigma_a_1;
nu_sigma_f_1=0.8;

thickness_region_1=1;

mesh_length_1=0.01;

%region 2

sigma_t_2=1;
sigma_a_2=0.5;
sigma_s_2=1-sigma_a_2;
nu_sigma_f_2=0;

thickness_region_2=1;

mesh_length_2=0.01;

%spatial discretization

x_region_1=(0:mesh_length_1:1)';
x_region_2=(1:mesh_length_2:2)';

edge_count_region_1=length(x_region_1);
mesh_count_region_1=edge_count_region_1-1;

edge_count_region_2=length(x_region_2);
mesh_count_region_2=edge_count_region_2-1;


x=cat(1,x_region_1(1:end-1,1),x_region_2(1:end,1));
edge_count=length(x);
mesh_count=edge_count-1;

mesh_length=zeros(mesh_count,1);
mesh_length(1:mesh_count_region_1)=mesh_length_1;
mesh_length(mesh_count_region_1+1:mesh_count)=mesh_length_2;


%data region specified data vector

vect_sigma_t=zeros(mesh_count,1);
vect_sigma_a=zeros(mesh_count,1);
vect_sigma_s=zeros(mesh_count,1);
vect_nu_sigma_f=zeros(mesh_count,1);

vect_sigma_t(1:mesh_count_region_1)=sigma_t_1;
vect_sigma_t(mesh_count_region_1+1:mesh_count)=sigma_t_2;


vect_sigma_a(1:mesh_count_region_1)=sigma_a_1;
vect_sigma_a(mesh_count_region_1+1:mesh_count)=sigma_a_2;

vect_sigma_s(1:mesh_count_region_1)=sigma_s_1;
vect_sigma_s(mesh_count_region_1+1:mesh_count)=sigma_s_2;

vect_nu_sigma_f(1:mesh_count_region_1)=nu_sigma_f_1;
vect_nu_sigma_f(mesh_count_region_1+1:mesh_count)=nu_sigma_f_2;

fission_term=vect_nu_sigma_f.*flux_old/(4*pi*k_old);
S=1/(4*pi)*vect_sigma_s.*flux_old+fission_term;
flux_new=transport_sweep(S);
iteration=1;

while(max(abs(flux_new-flux_old))>tol)
    flux_old=flux_new;
    S=1/(4*pi)*vect_sigma_s.*flux_old+fission_term;
    flux_new=transport_sweep(S);
    iteration=iteration+1;
end
iteration