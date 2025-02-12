function [potential,is_inside_solute,imaginary_part_of_potential]=calculate_potentials1(nodes,particles_params,particles_coefficients,medium_params,n_max,t_eps0_inv)
if size(nodes,2)~=3
    error('nodes must be an Nx3-array, where N is a number of nodes where the potential should be calculated.');    
end
potential=zeros([size(nodes,1),1]);
is_inside_solute=zeros([size(nodes,1),1]);
indices=1:size(nodes,1);
for i=1:numel(particles_params)
    ind=find(vecnorm(nodes-particles_params(i).center,2,2)<=particles_params(i).radius);
    if isempty(ind)==false
        is_inside_solute(indices(ind))=1;
        [phi,theta,tr]=cart2sph(nodes(ind,1)-particles_params(i).center(1),nodes(ind,2)-particles_params(i).center(2),nodes(ind,3)-particles_params(i).center(3));
        theta=pi/2-theta;
        coulombic_pot__excepting_center=particles_params(i).charge./(tr*particles_params(i).dielectric_constant);
        coulombic_pot__excepting_center(isinf(coulombic_pot__excepting_center))=0;
        potential(indices(ind))=calc_rn_Ynm(phi,cos(theta),medium_params.kappa*tr,particles_coefficients(i).L_nm,n_max)+coulombic_pot__excepting_center;
    end
    nodes(ind,:)=[];
    indices(ind)=[];
end
if isempty(indices)==false
    for i=1:numel(particles_params)
        [phi,theta,tr]=cart2sph(nodes(:,1)-particles_params(i).center(1),nodes(:,2)-particles_params(i).center(2),nodes(:,3)-particles_params(i).center(3));
        theta=pi/2-theta;
        tr=medium_params.kappa*tr;
        potential(indices)=potential(indices)+calc_kn_Ynm(phi,cos(theta),tr,particles_coefficients(i).G_nm,n_max);
    end
end

potential=t_eps0_inv*potential;
imaginary_part_of_potential=imag(potential);
potential=real(potential);
end

function res=calc_rn_Ynm(phi,cos_theta,tr,Li,n_max)
res=0;
incc=1;
temp_exp=exp(1i*(0:n_max).*phi);
for n=0:n_max
    temp=(tr.^n).*((((-1).^(0:n)).*(legendre(n,cos_theta,'norm')')/sqrt(2*pi)).*temp_exp(:,(0:n)+1));
    for m=(-n):(-1)
        res=res+Li(incc)*(conj(temp(:,1-m))*((-1)^(-m)));
        incc=incc+1;
    end
    for m=0:n
        res=res+Li(incc)*temp(:,1+m);
        incc=incc+1;
    end
end
end

function res=calc_kn_Ynm(phi,cos_theta,tr,Gi,n_max)
res=0;
incc=1;
temp_exp=exp(1i*(0:n_max).*phi);
for n=0:n_max
    temp=k_n(n,tr).*((((-1).^(0:n)).*(legendre(n,cos_theta,'norm')')/sqrt(2*pi)).*temp_exp(:,(0:n)+1));
    for m=(-n):(-1)
        res=res+Gi(incc)*(conj(temp(:,1-m))*((-1)^(-m)));
        incc=incc+1;
    end
    for m=0:n
        res=res+Gi(incc)*temp(:,1+m);
        incc=incc+1;
    end
end
end

function res=k_n(n,x)
res=sqrt(2/pi)*besselk(n+1/2,x)./sqrt(x);
end
