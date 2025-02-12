clear all
close all

medium_params.epsilon=80;

e_0  = 8.85418781762e-12; 
kb   = 1.380649e-23;    
T    = 273.15 + 25; 
e    = 1.602176634e-19;
N_av = 6.022e23;    
Angs = 1e-10;       
e_m  = 2.0;
e_w  = 80.;
ionic_strength = 0.145;
C_0 = 1.0e3 * N_av * ionic_strength; % Bulk concentration of monovalent species

eps_in  = 4.0*pi * e_0 * e_m * kb*T*Angs/(e^2); %adim e_in
eps_out = 4.0*pi * e_0 * e_w * kb*T*Angs/(e^2); %adim e_out

k = sqrt (2.0 * C_0 * Angs^2 * e^2 / (e_0 * e_w * kb * T));
l = 1/k;
t_eps0_inv = (e^2)/(4.0*pi * e_0 * kb*T*Angs);

medium_params.kappa= k;



fid = fopen ("../data/30spheres.pqr");
C = textscan (fid, '%*s %*s %*s %*s %*s %f %f %f %f %f', 'Delimiter', {' '}, 'MultipleDelimsAsOne', true, "CollectOutput", true);
pqr = C{1};
fclose (fid);
num_particles = size(pqr, 1);  % Number of rows
for i=1:num_particles
    particle_par.center=[pqr(i,1),pqr(i,2),pqr(i,3)];
    particle_par.charge=pqr(i,4);
    particle_par.radius=pqr(i,5);
    particle_par.dielectric_constant=2;    
    particles_params(i)=particle_par;
end
%YOU MAY ADD HERE (IN THE SIMILAR MANNER) ARBITRARY NUMBER OF OTHER SPHERES

load array_cg15.mat array_clebschgordan cg_n_max;

n_max=15;
tic;
%calculate the total energy
tol_gmres=1e-12; %tolerance for the gmres iterative solver (to determine potential expansion coefficients "particles_coefficients", see below)
[particles_coefficients,energy]=multi_spheres_cg5(particles_params,medium_params,n_max,array_clebschgordan,cg_n_max,tol_gmres,t_eps0_inv);
toc;
tic;
%calculate the energy components
[energy_reaction,energy_ionic,energy_Coulombic] = calc_energy_components0(particles_params,medium_params,n_max,particles_coefficients,energy,t_eps0_inv);
toc;
tic;
%calculate the potential specic point
%positions_pot =load ("path/to/your/data/points");
xmin=-20;
xmax=20;
ymin=-20;
ymax=20;
zmin=-20;
zmax=70;
step=2;
[y,z,x]=meshgrid(xmin:step:xmax,zmin:step:zmax,ymin:step:ymax);
nodes=[x(:),y(:),z(:)];
[potential,is_inside_solute2,dummy_image_part2]=calculate_potentials1(nodes,particles_params,particles_coefficients,medium_params,n_max,t_eps0_inv);
toc;

[Xid,Yid,Zid] = sphere;
for i=1:numel(particles_params)
    X=Xid*particles_params(i).radius;
    Y=Yid*particles_params(i).radius;
    Z=Zid*particles_params(i).radius;
    surf(X+particles_params(i).center(1),Y+particles_params(i).center(2),Z+particles_params(i).center(3));
    hold on;
end
axis equal;
xlabel('X');ylabel('Y');zlabel('Z');