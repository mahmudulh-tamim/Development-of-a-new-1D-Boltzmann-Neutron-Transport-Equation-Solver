function flux=transport_sweep(S,group_i)




%% group data
sigma_t=zeros(6,1);
sigma_f=zeros(6,1);
chi=zeros(6,1);
nu=zeros(6,1);

sigma_s=zeros(6,6);

%group 1

sigma_t(1,1)=0.24;
sigma_f(1,1)=0.006;
chi(1,1)=0.48;
nu(1,1)=3;
sigma_s(1,1:6)=[0.024, 0.0, 0.0, 0.0, 0.0, 0.0];

%group 2

chi(2,1)=0.02;
nu(2,1)=2.5;
sigma_f(2,1)=0.06;
sigma_t(2,1)=0.975;
sigma_s(2,1:6)=[0.171, 0.6, 0.0, 0.0, 0.0, 0.0];

%group 3

chi(3,1)=0.0;
nu(3,1)=2;
sigma_f(3,1)=0.9;
sigma_t(3,1)=3.10;
sigma_s(3,1:6)=[0.033, 0.275, 2.0, 0.0, 0.0, 0.0];

%group 4

chi(4,1)=0.0;
nu(4,1)=2;
sigma_f(4,1)=0.9;
sigma_t(4,1)=3.10;
sigma_s(4,1:6)=[0.0, 0.0, 0.0, 2.0, 0.275, 0.033];

%group 5

chi(5,1)=0.02;
nu(5,1)=2.5;
sigma_f(5,1)=0.06;
sigma_t(5,1)=0.975;
sigma_s(5,1:6)=[0.0, 0.0, 0.0, 0.0, 0.60, 0.171];

%group 6

sigma_t(6,1)=0.24;
sigma_f(6,1)=0.006;
chi(6,1)=0.48;
nu(6,1)=3;
sigma_s(6,1:6)=[0.0, 0.0, 0.0, 0.0, 0.0, 0.024];


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







%% initialization

%initialization of angular flux

psi_out=zeros(polar_discretization_number, azimuthal_discretization_number,edge_count);

del_psi=zeros(polar_discretization_number, azimuthal_discretization_number,mesh_count);

avg_psi=zeros(polar_discretization_number, azimuthal_discretization_number,mesh_count);

psi_out_old=psi_out;
%ray tracing length

traced_ray=zeros(polar_discretization_number, azimuthal_discretization_number,mesh_count);

%scaler flux

flux=zeros(mesh_count,1);

%ray tracing left to right

for i=1:edge_count
    for j=1:azimuthal_discretization_number/2
        for p=1:polar_discretization_number
            
            if(i==1)
                psi_out(p,j,i)=0;
            else
                traced_ray(p,j,i-1)=mesh_length(i-1,1)/abs(mu(p,1)*cos(azimuthal_direction_theta(j,1)));

                del_psi(p,j,i-1)=(psi_out(p,j,i-1)-S(i-1,1)/sigma_t(group_i,1))*(1-exp(-sigma_t(group_i,1)*traced_ray(p,j,i-1)));
                
                psi_out(p,j,i)=psi_out(p,j,i-1)-del_psi(p,j,i-1);

                avg_psi(p,j,i-1)=del_psi(p,j,i-1)/(traced_ray(p,j,i-1)*sigma_t(group_i,1))+S(i-1,1)/sigma_t(group_i,1);
               % flux(i-1,1)=flux(i-1,1)+avg_psi(p,j,i-1)*del_theta*w(p,1);
            end
        end
    end
end

%ray tracing right to left

for i=edge_count:-1:1
    for j=azimuthal_discretization_number/2+1:azimuthal_discretization_number
        for p=1:polar_discretization_number
            if(i==edge_count)
                psi_out(p,j,i)=psi_out(p,j-azimuthal_discretization_number/2,i);
            else
                traced_ray(p,j,i)=mesh_length(i,1)/abs(mu(p,1)*cos(azimuthal_direction_theta(j,1)));
                del_psi(p,j,i)=(psi_out(p,j,i+1)-S(i,1)/sigma_t(group_i,1))*(1-exp(-sigma_t(group_i,1)*traced_ray(p,j,i)));
                psi_out(p,j,i)=psi_out(p,j,i+1)-del_psi(p,j,i);

                avg_psi(p,j,i)=del_psi(p,j,i)/(traced_ray(p,j,i)*sigma_t(group_i,1))+S(i,1)/sigma_t(group_i,1);
               % flux(i,1)=flux(i,1)+avg_psi(p,j,i)*del_theta*w(p,1);
            end
            
                
        end
    end
end

transport_iteration=1;
ghojad=reshape(abs(psi_out-psi_out_old),[],1);
while(max(ghojad)>10^(-6))
    psi_out_old=psi_out;
        %ray tracing left to right

        for i=1:edge_count
            for j=1:azimuthal_discretization_number/2
                for p=1:polar_discretization_number
                    
                    if(i==1)
                        psi_out(p,j,i)=psi_out(p,j+azimuthal_discretization_number/2,i);
                    else
                        traced_ray(p,j,i-1)=mesh_length(i-1,1)/abs(mu(p,1)*cos(azimuthal_direction_theta(j,1)));
        
                        del_psi(p,j,i-1)=(psi_out(p,j,i-1)-S(i-1,1)/sigma_t(group_i,1))*(1-exp(-sigma_t(group_i,1)*traced_ray(p,j,i-1)));
                        
                        psi_out(p,j,i)=psi_out(p,j,i-1)-del_psi(p,j,i-1);
        
                        avg_psi(p,j,i-1)=del_psi(p,j,i-1)/(traced_ray(p,j,i-1)*sigma_t(group_i,1))+S(i-1,1)/sigma_t(group_i,1);
                       % flux(i-1,1)=flux(i-1,1)+avg_psi(p,j,i-1)*del_theta*w(p,1);
                    end
                end
            end
        end
        
        %ray tracing right to left
        
        for i=edge_count:-1:1
            for j=azimuthal_discretization_number/2+1:azimuthal_discretization_number
                for p=1:polar_discretization_number
                    if(i==edge_count)
                        psi_out(p,j,i)=psi_out(p,j-azimuthal_discretization_number/2,i);
                    else
                        traced_ray(p,j,i)=mesh_length(i,1)/abs(mu(p,1)*cos(azimuthal_direction_theta(j,1)));
                        del_psi(p,j,i)=(psi_out(p,j,i+1)-S(i,1)/sigma_t(group_i,1))*(1-exp(-sigma_t(group_i,1)*traced_ray(p,j,i)));
                        psi_out(p,j,i)=psi_out(p,j,i+1)-del_psi(p,j,i);
        
                        avg_psi(p,j,i)=del_psi(p,j,i)/(traced_ray(p,j,i)*sigma_t(group_i,1))+S(i,1)/sigma_t(group_i,1);
                       % flux(i,1)=flux(i,1)+avg_psi(p,j,i)*del_theta*w(p,1);
                    end
                    
                        
                end
            end
        end
        transport_iteration=1+transport_iteration;
        ghojad=reshape(abs(psi_out-psi_out_old),[],1);

end

%ray tracing left to right

for i=2:edge_count
    for j=1:azimuthal_discretization_number/2
        for p=1:polar_discretization_number
              flux(i-1,1)=flux(i-1,1)+avg_psi(p,j,i-1)*del_theta*w(p,1);
           
        end
    end
end

%ray tracing right to left

for i=edge_count-1:-1:1
    for j=azimuthal_discretization_number/2+1:azimuthal_discretization_number
        for p=1:polar_discretization_number
            
            flux(i,1)=flux(i,1)+avg_psi(p,j,i)*del_theta*w(p,1);
         end
          
     end
 end

transport_iteration
