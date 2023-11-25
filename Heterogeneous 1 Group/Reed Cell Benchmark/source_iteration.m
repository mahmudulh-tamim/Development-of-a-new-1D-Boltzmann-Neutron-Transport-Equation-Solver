function source_iteration(flux_old)

tol=10^(-6);

%region 1
sigma_t_1=50;
sigma_s_1=0;
S_1=40;


thickness_1=2;
mesh_length_1=thickness_1/200;

%region 2

sigma_t_2=5;
sigma_s_2=0;
S_2=0;

thickness_2=1;
mesh_length_2=thickness_2/100;

%region 3
sigma_t_3=0;
sigma_s_3=0;
S_3=0;

thickness_3=2;
mesh_length_3=thickness_3/200;

%region 4
sigma_t_4=1;
sigma_s_4=0.9;
S_4=1;

thickness_4=1;
mesh_length_4=thickness_4/100;

%region 5
sigma_t_5=1;
sigma_s_5=0.9;
S_5=0;

thickness_5=2;
mesh_length_5=thickness_5/200;




%spatial discretization

mesh_length_1=mesh_length_1;
mesh_length_2=mesh_length_2;
mesh_length_3=mesh_length_3;
mesh_length_4=mesh_length_4;
mesh_length_5=mesh_length_5;


x_region_1=(0:mesh_length_1:thickness_1)';
x_region_2=(thickness_1:mesh_length_2:thickness_1+thickness_2)';
x_region_3=(thickness_1+thickness_2:mesh_length_3:thickness_1+thickness_2+thickness_3)';
x_region_4=(thickness_1+thickness_2+thickness_3:mesh_length_4:thickness_1+thickness_2+thickness_3+thickness_4)';
x_region_5=(thickness_1+thickness_2+thickness_3+thickness_4:mesh_length_5:thickness_1+thickness_2+thickness_3+thickness_4+thickness_5)';




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




x=cat(1,x_region_1(1:end-1,1),x_region_2(1:end-1,1),x_region_3(1:end-1,1), x_region_4(1:end-1,1),x_region_5(1:end,1));
edge_count=length(x);
mesh_count=edge_count-1;

mesh_length=zeros(mesh_count,1);

mesh_length(1:mesh_count_region_1,1)=mesh_length_1;
mesh_length(mesh_count_region_1+1:mesh_count_region_2,1)=mesh_length_2;
mesh_length(mesh_count_region_2+1:mesh_count_region_3,1)=mesh_length_3;
mesh_length(mesh_count_region_3+1:mesh_count_region_4,1)=mesh_length_4;
mesh_length(mesh_count_region_4+1:mesh_count_region_5,1)=mesh_length_5;


%region specified data vector
vect_sigma_t=zeros(mesh_count,1);
vect_sigma_s=zeros(mesh_count,1);
Ss=zeros(mesh_count,1);

vect_sigma_t(1:mesh_count_region_1,1)=sigma_t_1;
vect_sigma_t(mesh_count_region_1+1:mesh_count_region_2,1)=sigma_t_2;
vect_sigma_t(mesh_count_region_2+1:mesh_count_region_3,1)=sigma_t_3;
vect_sigma_t(mesh_count_region_3+1:mesh_count_region_4,1)=sigma_t_4;
vect_sigma_t(mesh_count_region_4+1:mesh_count_region_5,1)=sigma_t_5;


vect_sigma_s(1:mesh_count_region_1,1)=sigma_s_1;
vect_sigma_s(mesh_count_region_1+1:mesh_count_region_2,1)=sigma_s_2;
vect_sigma_s(mesh_count_region_2+1:mesh_count_region_3,1)=sigma_s_3;
vect_sigma_s(mesh_count_region_3+1:mesh_count_region_4,1)=sigma_s_4;
vect_sigma_s(mesh_count_region_4+1:mesh_count_region_5,1)=sigma_s_5;


Ss(1:mesh_count_region_1,1)=S_1;
Ss(mesh_count_region_1+1:mesh_count_region_2,1)=S_2;
Ss(mesh_count_region_2+1:mesh_count_region_3,1)=S_3;
Ss(mesh_count_region_3+1:mesh_count_region_4,1)=S_4;
Ss(mesh_count_region_4+1:mesh_count_region_5,1)=S_5;





S=1/(4*pi)*vect_sigma_s.*flux_old+Ss/(4*pi);
flux_new=transport_sweep(S);
iteration=1;

while(max(abs(flux_new-flux_old))>tol)
    flux_old=flux_new;
    S=1/(4*pi)*vect_sigma_s.*flux_old+Ss/(4*pi);
    flux_new=transport_sweep(S);
    iteration=iteration+1;
end
iteration

x_mid=0.5*(x(1:end-1)+x(2:end));

plot(x_mid, flux_new);

max(flux_new)
