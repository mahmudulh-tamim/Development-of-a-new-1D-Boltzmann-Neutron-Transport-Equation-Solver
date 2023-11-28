function flux_new_ith_group=source_iteration(flux_old,k_old,group_i)

tol=10^(-6);

%initialization
chi_U=zeros(1,2);
nu_U=zeros(1,2);
sigma_f_U=zeros(1,2);
sigma_t_U=zeros(1,2);
sigma_s_U=zeros(1,2,2);

chi_H2O=zeros(1,2);
nu_H2O=zeros(1,2);
sigma_f_H2O=zeros(1,2);
sigma_t_H2O=zeros(1,2);
sigma_s_H2O=zeros(1,2,2);

%% group data
%uranium
%group 1

chi_U(1,1)=1;
nu_U(1,1)=1.004;
sigma_f_U(1,1)=0.61475;
sigma_t_U(1,1)=0.650917;
sigma_s_U(1,:,1)=[0,0];

%group 2

chi_U(1,2)=0;
nu_U(1,2)=2.5;
sigma_f_U(1,2)=0.045704;
sigma_t_U(1,2)=2.138;
sigma_s_U(1,:,2)=[0.0342008,2.0688];

%H2O
%group 1

chi_H2O(1,1)=0;
nu_H2O(1,1)=0;
sigma_f_H2O(1,1)=0;
sigma_t_H2O(1,1)=1.331518007;
sigma_s_H2O(1,:,1)=[1.226381244,0];

%group 2

chi_H2O(1,2)=0;
nu_H2O(1,2)=0;
sigma_f_H2O(1,2)=0;
sigma_t_H2O(1,2)=4.3735;
sigma_s_H2O(1,:,2)=[0.104639534,4.3547];

%% geometry data
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


%% angular discretization

%polar discretization

mu=[0.932954;0.537707;0.166648;-0.166648;-0.537707;-0.932954];
w=[0.670148;0.283619;0.046233;0.046233;0.283619;0.670148];
polar_discretization_number=size(mu,1);

%azimuthal discretization
N_a=128;
del_theta=2*pi/N_a;
theta=(0:del_theta:2*pi)';
azimuthal_direction_theta= 0.5*(theta(1:end-1,1)+theta(2:end,1));
azimuthal_discretization_number=size(azimuthal_direction_theta,1);

%end of geometry data


%% region and group specified  material properties

chi=zeros(mesh_count,2);
nu=zeros(mesh_count,2);
sigma_f=zeros(mesh_count,2);
sigma_t=zeros(mesh_count,2);

sigma_s=zeros(mesh_count,2,2);

%chi
chi(1:half_mesh_count_H2O,1)=chi_H2O(1,1);
chi(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,1)=chi_U(1,1);
chi(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,1)=chi_U(1,1);
chi(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,1)=chi_H2O(1,1);

chi(1:half_mesh_count_H2O,2)=chi_H2O(1,2);
chi(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,2)=chi_U(1,2);
chi(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,2)=chi_U(1,2);
chi(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,2)=chi_H2O(1,2);

%nu
nu(1:half_mesh_count_H2O,1)=nu_H2O(1,1);
nu(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,1)=nu_U(1,1);
nu(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,1)=nu_U(1,1);
nu(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,1)=nu_H2O(1,1);


nu(1:half_mesh_count_H2O,2)=nu_H2O(1,2);
nu(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,2)=nu_U(1,2);
nu(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,2)=nu_U(1,2);
nu(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,2)=nu_H2O(1,2);

%sigma_f
sigma_f(1:half_mesh_count_H2O,1)=sigma_f_H2O(1,1);
sigma_f(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,1)=sigma_f_U(1,1);
sigma_f(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,1)=sigma_f_U(1,1);
sigma_f(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,1)=sigma_f_H2O(1,1);


sigma_f(1:half_mesh_count_H2O,2)=sigma_f_H2O(1,2);
sigma_f(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,2)=sigma_f_U(1,2);
sigma_f(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,2)=sigma_f_U(1,2);
sigma_f(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,2)=sigma_f_H2O(1,2);

%sigma_t
sigma_t(1:half_mesh_count_H2O,1)=sigma_t_H2O(1,1);
sigma_t(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,1)=sigma_t_U(1,1);
sigma_t(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,1)=sigma_t_U(1,1);
sigma_t(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,1)=sigma_t_H2O(1,1);


sigma_t(1:half_mesh_count_H2O,2)=sigma_t_H2O(1,2);
sigma_t(half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U,2)=sigma_t_U(1,2);
sigma_t(half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U,2)=sigma_t_U(1,2);
sigma_t(half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U,2)=sigma_t_H2O(1,2);


%sigma_s
for h=1:half_mesh_count_H2O
    sigma_s(h,:,:)=sigma_s_H2O;
end
for h=half_mesh_count_H2O+1:half_mesh_count_H2O+half_mesh_count_U
    sigma_s(h,:,:)=sigma_s_U;
end
for h=half_mesh_count_H2O+half_mesh_count_U+1:half_mesh_count_H2O+2*half_mesh_count_U
    sigma_s(h,:,:)=sigma_s_U;
end
for h=half_mesh_count_H2O+2*half_mesh_count_U+1:2*half_mesh_count_H2O+2*half_mesh_count_U
    sigma_s(h,:,:)=sigma_s_H2O;
end


%% source iteration

flux_new=flux_old;

fission_term=chi(:,group_i).*(nu.*sigma_f).*flux_old/(4*pi*k_old);
real_fission_term=sum(fission_term,2);

scattering_term=sigma_s(:,:,group_i).*flux_old/(4*pi);
real_scattering_term=sum(scattering_term,2);


S=real_fission_term+real_scattering_term;
flux_new(:,group_i)=transport_sweep(S,group_i);

iteration=1;

while(max(abs(flux_new(:,group_i)-flux_old(:,group_i)))>tol)

    flux_old(:,group_i)=flux_new(:,group_i);
    scattering_term=sigma_s(:,:,group_i).*flux_old/(4*pi);
    real_scattering_term=sum(scattering_term,2);
    S=real_fission_term+real_scattering_term;
    flux_new(:,group_i)=transport_sweep(S,group_i);
    iteration=1+iteration;
end
iteration
flux_new_ith_group=flux_new(:,group_i);









