function [energy_reaction,energy_ionic,energy_Coulombic] = calc_energy_components0(particles_params,medium_params,n_max,particles_coefficients,energy_total, t_eps0_inv)
% the parameter "energy_total" (the total electrostatic energy) is calculated from "multi_spheres_cg...m"
kappa=medium_params.kappa;
epsilon_m=medium_params.epsilon;
%REACTION ENERGY:
energy_reaction=0;
for i=1:numel(particles_params)
    tai=kappa*particles_params(i).radius;
    epsilon_i=particles_params(i).dielectric_constant;
    energy_reaction=energy_reaction+particles_params(i).charge*particles_params(i).charge*(1/epsilon_m-1/epsilon_i)/particles_params(i).radius;
    for j=[1:(i-1),(i+1):numel(particles_params)]
%         Rij_vec=particles_params(j).center-particles_params(i).center;
        [phi,theta,Rij]=cart2sph(particles_params(j).center(1)-particles_params(i).center(1),particles_params(j).center(2)-particles_params(i).center(2),particles_params(j).center(3)-particles_params(i).center(3));
        theta=pi/2-theta;
%         Rij=sqrt(Rij_vec(:,1).^2+Rij_vec(:,2).^2+Rij_vec(:,3).^2);
        Phi_polarization_ij=particles_params(i).charge/Rij;
        temp_exp=exp(1i*(0:n_max).*phi);
        incc=2;%because we have factor $n$, so we can start summation from n=1
        for n=1:n_max
            temp_Ynm=(((-1).^(0:n)).*(legendre(n,cos(theta),'norm')')/sqrt(2*pi)).*temp_exp(:,(0:n)+1);
            temp_m=0;
            for m=(-n):(-1)
                temp_m=temp_m+particles_coefficients(i).L_nm(incc)*(conj(temp_Ynm(:,1-m))*((-1)^(-m)));
                incc=incc+1;
            end
            for m=0:n
                temp_m=temp_m+particles_coefficients(i).L_nm(incc)*temp_Ynm(:,1+m);
                incc=incc+1;
            end
            Phi_polarization_ij=Phi_polarization_ij-epsilon_i*(n/(2*n+1))*(tai^n)*((particles_params(i).radius/Rij)^(n+1))*temp_m;
        end
        Phi_polarization_ij=Phi_polarization_ij*(1/epsilon_m-1/epsilon_i);
        energy_reaction=energy_reaction+particles_params(j).charge*Phi_polarization_ij;
    end
end


energy_reaction=(t_eps0_inv*energy_reaction/2);
%COULOMBIC ENERGY:
energy_Coulombic=0;
for i=1:numel(particles_params)
    temp_pot=0;
    for j=[1:(i-1),(i+1):numel(particles_params)]
        Rij_vec=particles_params(j).center-particles_params(i).center;
        temp_pot=temp_pot+particles_params(j).charge/(particles_params(j).dielectric_constant*sqrt(Rij_vec(:,1).^2+Rij_vec(:,2).^2+Rij_vec(:,3).^2));
    end
    energy_Coulombic=energy_Coulombic+particles_params(i).charge*temp_pot;
end
energy_Coulombic=(t_eps0_inv*energy_Coulombic/2);
energy_ionic=(energy_total-energy_reaction-energy_Coulombic);
end