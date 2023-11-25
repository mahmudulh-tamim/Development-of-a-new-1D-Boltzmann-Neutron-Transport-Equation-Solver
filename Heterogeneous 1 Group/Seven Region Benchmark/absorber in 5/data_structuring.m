%reflector
sigma_t_reflector=0.371;
sigma_s_reflector=0.334;
nu_sigma_f_reflector=0;

thickness_reflector=1/sigma_t_reflector;
mesh_length_reflector=thickness_reflector/125;

%fuel

sigma_t_fuel=0.415;
sigma_s_fuel=0.334;
nu_sigma_f_fuel=0.178;

thickness_fuel=1/sigma_t_fuel;
mesh_length_fuel=thickness_fuel/125;

%absorber
sigma_t_absorber=0.371;
sigma_s_absorber=0.037;
nu_sigma_f_absorber=0;

thickness_absorber=1/sigma_t_absorber;
mesh_length_absorber=thickness_absorber/125;


%spatial discretization

mesh_length_1=mesh_length_reflector;
mesh_length_2=mesh_length_fuel;
mesh_length_3=mesh_length_reflector;
mesh_length_4=mesh_length_fuel;
mesh_length_5=mesh_length_absorber;
mesh_length_6=mesh_length_fuel;
mesh_length_7=mesh_length_reflector;

x_region_1=(0:mesh_length_reflector:1*thickness_reflector)';
x_region_2=(1*thickness_reflector:mesh_length_fuel:1*thickness_reflector+1*thickness_fuel)';
x_region_3=(1*thickness_reflector+1*thickness_fuel:mesh_length_reflector:2*thickness_reflector+1*thickness_fuel)';
x_region_4=(2*thickness_reflector+1*thickness_fuel:mesh_length_fuel:2*thickness_reflector+2*thickness_fuel)';
x_region_5=(2*thickness_reflector+2*thickness_fuel:mesh_length_absorber:2*thickness_reflector+2*thickness_fuel+thickness_absorber)';
x_region_6=(2*thickness_reflector+2*thickness_fuel+thickness_absorber:mesh_length_fuel:2*thickness_reflector+3*thickness_fuel+thickness_absorber)';
x_region_7=(2*thickness_reflector+3*thickness_fuel+thickness_absorber:mesh_length_reflector:3*thickness_reflector+3*thickness_fuel+thickness_absorber)';




edge_count_region_1=length(x_region_1);
mesh_count_region_1=edge_count_region_1-1;

edge_count_region_2=length(x_region_2);
mesh_count_region_2=edge_count_region_2-1+mesh_count_region_1;

edge_count_region_3=length(x_region_3);
mesh_count_region_3=edge_count_region_3-1+mesh_count_region_2;

edge_count_region_4=length(x_region_4);
mesh_count_region_4=edge_count_region_4-1+mesh_count_region_3;

edge_count_region_5=length(x_region_5);
mesh_count_region_5=edge_count_region_5-1+mesh_count_region_4;

edge_count_region_6=length(x_region_6);
mesh_count_region_6=edge_count_region_6-1+mesh_count_region_5;

edge_count_region_7=length(x_region_7);
mesh_count_region_7=edge_count_region_7-1+mesh_count_region_6;




x=cat(1,x_region_1(1:end-1,1),x_region_2(1:end-1,1),x_region_3(1:end-1,1), x_region_4(1:end-1,1),x_region_5(1:end-1,1),x_region_6(1:end-1,1),x_region_7(1:end,1));
edge_count=length(x);
mesh_count=edge_count-1;

mesh_length=zeros(mesh_count,1);

mesh_length(1:mesh_count_region_1,1)=mesh_length_1;
mesh_length(mesh_count_region_1+1:mesh_count_region_2,1)=mesh_length_2;
mesh_length(mesh_count_region_2+1:mesh_count_region_3,1)=mesh_length_3;
mesh_length(mesh_count_region_3+1:mesh_count_region_4,1)=mesh_length_4;
mesh_length(mesh_count_region_4+1:mesh_count_region_5,1)=mesh_length_5;
mesh_length(mesh_count_region_5+1:mesh_count_region_6,1)=mesh_length_6;
mesh_length(mesh_count_region_6+1:mesh_count_region_7,1)=mesh_length_7;

%region specified data vector
vect_sigma_t=zeros(mesh_count,1);
vect_sigma_s=zeros(mesh_count,1);
vect_nu_sigma_f=zeros(mesh_count,1);

vect_sigma_t(1:mesh_count_region_1,1)=sigma_t_reflector;
vect_sigma_t(mesh_count_region_1+1:mesh_count_region_2,1)=sigma_t_fuel;
vect_sigma_t(mesh_count_region_2+1:mesh_count_region_3,1)=sigma_t_reflector;
vect_sigma_t(mesh_count_region_3+1:mesh_count_region_4,1)=sigma_t_fuel;
vect_sigma_t(mesh_count_region_4+1:mesh_count_region_5,1)=sigma_t_absorber;
vect_sigma_t(mesh_count_region_5+1:mesh_count_region_6,1)=sigma_t_fuel;
vect_sigma_t(mesh_count_region_6+1:mesh_count_region_7,1)=sigma_t_reflector;

vect_sigma_s(1:mesh_count_region_1,1)=sigma_s_reflector;
vect_sigma_s(mesh_count_region_1+1:mesh_count_region_2,1)=sigma_s_fuel;
vect_sigma_s(mesh_count_region_2+1:mesh_count_region_3,1)=sigma_s_reflector;
vect_sigma_s(mesh_count_region_3+1:mesh_count_region_4,1)=sigma_s_fuel;
vect_sigma_s(mesh_count_region_4+1:mesh_count_region_5,1)=sigma_s_absorber;
vect_sigma_s(mesh_count_region_5+1:mesh_count_region_6,1)=sigma_s_fuel;
vect_sigma_s(mesh_count_region_6+1:mesh_count_region_7,1)=sigma_s_reflector;

vect_nu_sigma_f(1:mesh_count_region_1,1)=nu_sigma_f_reflector;
vect_nu_sigma_f(mesh_count_region_1+1:mesh_count_region_2,1)=nu_sigma_f_fuel;
vect_nu_sigma_f(mesh_count_region_2+1:mesh_count_region_3,1)=nu_sigma_f_reflector;
vect_nu_sigma_f(mesh_count_region_3+1:mesh_count_region_4,1)=nu_sigma_f_fuel;
vect_nu_sigma_f(mesh_count_region_4+1:mesh_count_region_5,1)=nu_sigma_f_absorber;
vect_nu_sigma_f(mesh_count_region_5+1:mesh_count_region_6,1)=nu_sigma_f_fuel;
vect_nu_sigma_f(mesh_count_region_6+1:mesh_count_region_7,1)=nu_sigma_f_reflector;

%% angular discretization

%polar discretization

mu=[0.932954;0.537707;0.166648;-0.166648;-0.537707;-0.932954];
w=[0.670148;0.283619;0.046233;0.046233;0.283619;0.670148];
polar_discretization_number=size(mu,1);

%azimuthal discretization
N_a=128;
del_nu=2*pi/N_a;
nu=(0:del_nu:2*pi)';
azimuthal_direction_nu= 0.5*(nu(1:end-1,1)+nu(2:end,1));
azimuthal_discretization_number=size(azimuthal_direction_nu,1);
