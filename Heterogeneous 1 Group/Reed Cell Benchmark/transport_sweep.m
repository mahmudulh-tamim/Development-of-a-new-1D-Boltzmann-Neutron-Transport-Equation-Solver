

function flux=transport_sweep(S)

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


%% initialization

%initialization of angular flux

psi_out=zeros(polar_discretization_number, azimuthal_discretization_number,edge_count);

del_psi=zeros(polar_discretization_number, azimuthal_discretization_number,mesh_count);

avg_psi=zeros(polar_discretization_number, azimuthal_discretization_number,mesh_count);

%ray tracing length

traced_ray=zeros(polar_discretization_number, azimuthal_discretization_number,mesh_count);

%scaler flux

flux=zeros(mesh_count,1);
%ray tracing right to left

for i=edge_count:-1:1
    for j=azimuthal_discretization_number/2+1:azimuthal_discretization_number
        for p=1:polar_discretization_number
            if(i<301||i>501)
                if(i==edge_count)
                    psi_out(p,j,i)=0;
                else
                    traced_ray(p,j,i)=mesh_length(i)/abs(mu(p,1)*cos(azimuthal_direction_nu(j,1)));
                    del_psi(p,j,i)=(psi_out(p,j,i+1)-S(i,1)/vect_sigma_t(i,1))*(1-exp(-vect_sigma_t(i,1)*traced_ray(p,j,i)));
                    psi_out(p,j,i)=psi_out(p,j,i+1)-del_psi(p,j,i);
    
                    avg_psi(p,j,i)=del_psi(p,j,i)/(traced_ray(p,j,i)*vect_sigma_t(i,1))+S(i,1)/vect_sigma_t(i,1);
                    flux(i,1)=flux(i,1)+avg_psi(p,j,i)*del_nu*w(p,1);
                end
            else
                traced_ray(p,j,i)=mesh_length(i)/abs(mu(p,1)*cos(azimuthal_direction_nu(j,1)));
                del_psi(p,j,i)=(psi_out(p,j,i+1))*(1-exp(-vect_sigma_t(i,1)*traced_ray(p,j,i)));
                psi_out(p,j,i)=psi_out(p,j,i+1)-del_psi(p,j,i);

                avg_psi(p,j,i)=psi_out(p,j,i);
                flux(i,1)=flux(i,1)+avg_psi(p,j,i)*del_nu*w(p,1);
            end
        end
    end
end

%ray tracing left to right

for i=1:edge_count
    for j=1:azimuthal_discretization_number/2
        for p=1:polar_discretization_number
            if(i<=301||i>=502)
                if(i==1)
                    psi_out(p,j,i)=psi_out(p,j+azimuthal_discretization_number/2,i);
                else
                    traced_ray(p,j,i-1)=mesh_length(i-1)/abs(mu(p,1)*cos(azimuthal_direction_nu(j,1)));
    
                    del_psi(p,j,i-1)=(psi_out(p,j,i-1)-S(i-1,1)/vect_sigma_t(i-1,1))*(1-exp(-vect_sigma_t(i-1,1)*traced_ray(p,j,i-1)));
                    
                    psi_out(p,j,i)=psi_out(p,j,i-1)-del_psi(p,j,i-1);
    
                    avg_psi(p,j,i-1)=del_psi(p,j,i-1)/(traced_ray(p,j,i-1)*vect_sigma_t(i-1,1))+S(i-1,1)/vect_sigma_t(i-1,1);
                    flux(i-1,1)=flux(i-1,1)+avg_psi(p,j,i-1)*del_nu*w(p,1);
                end
            else
                 traced_ray(p,j,i-1)=mesh_length(i-1)/abs(mu(p,1)*cos(azimuthal_direction_nu(j,1)));
    
                 del_psi(p,j,i-1)=(psi_out(p,j,i-1))*(1-exp(-vect_sigma_t(i-1,1)*traced_ray(p,j,i-1)));
                    
                 psi_out(p,j,i)=psi_out(p,j,i-1)-del_psi(p,j,i-1);
    
                 avg_psi(p,j,i-1)=psi_out(p,j,i);
                 flux(i-1,1)=flux(i-1,1)+avg_psi(p,j,i-1)*del_nu*w(p,1);
            end

        end
    end
end


end

