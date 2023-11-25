function power_iteration()
tol=10^(-12);

%remember nu has a changed meaning now so change the variable name for
%polar angle

%% group data
sigma_t=zeros(3,1);
sigma_f=zeros(3,1);
chi=zeros(3,1);
nu=zeros(3,1);

sigma_s=zeros(3,3);

%group 1

sigma_t(1,1)=0.24;
sigma_f(1,1)=0.006;
chi(1,1)=0.96;
nu(1,1)=3;
sigma_s(1,1:3)=[0.024, 0.0, 0.0];

%group 2

sigma_t(2,1)=0.975;
sigma_f(2,1)=0.06;
chi(2,1)=0.04;
nu(2,1)=2.5;
sigma_s(2,1:3)=[0.171, 0.6, 0.0];

%group 3

sigma_t(3,1)=3.1;
sigma_f(3,1)=0.9;
chi(3,1)=0.0;
nu(3,1)=2;
sigma_s(3,1:3)=[0.033, 0.275, 2.0];

%% end of group data

%% geometry data
%spatial discretization

slab_length=2;
mesh_number=100;
del_x=2/mesh_number;
x=(0:del_x:2)';
edge_count=length(x);
mesh_count=edge_count-1; % it should be equal to mesh_number

mesh_length=zeros(mesh_count,1);
mesh_length(1:end,1)=del_x;

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

%% end of geometry data





k_old=1;
flux_old=ones(mesh_count,3);%/sum((mesh_length.*vect_nu_sigma_f))

flux_new_iteration=group_flux(flux_old,k_old);

k_new=k_old*sum((nu.*sigma_f)'.*(mesh_length.*flux_new_iteration))/sum((nu.*sigma_f)'.*(mesh_length.*flux_old));



power_iteration_count=1;

while (max(abs(flux_new_iteration-flux_old))>tol)
    flux_old=flux_new_iteration;
    k_old=k_new;
    flux_new_iteration=group_flux(flux_old,k_old);

    k_new=k_old*sum((nu.*sigma_f)'.*(mesh_length.*flux_new_iteration))/sum((nu.*sigma_f)'.*(mesh_length.*flux_old));

    power_iteration_count=power_iteration_count+1;
end

x_mid=0.5*(x(1:end-1)+x(2:end));

plot(x_mid, flux_new_iteration);

power_iteration_count
k_new
flux_new_iteration(:,2)./flux_new_iteration(:,1)

flux_new_iteration(:,3)./flux_new_iteration(:,1)
flux_new_iteration(:,3)./flux_new_iteration(:,2)
