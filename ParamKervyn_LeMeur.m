%% Matlab
clear; close all; clc;
Nsteps = 1000;

% Parametres physiques:
tfin = 20.0;		         % temps physique de la simulation (s)
xL   = 0.0 ; 		     % limite gauche du fond marin (m) const: 0.
xR   = 10.0 ; 		     % limite droite du fond marin (m) const: 15.
yL   = 0.0 ;		     % limite inferieur du fond marin (m) const: 0.
yU   = 6.0 ;		     % limite superieur du fond marin (m) const: 15.
pert_amplitude = 0.e0;	 % amplitude de la perturbation du terme de droite (m)
pert_velocity  = 0.e0;	 % vitesse angulaire de la perturbation du terme de droite (m)

% Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

% Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 2.e0;      % profondeur maximale du fond marin (m)
h1  = 1.e0;      % profondeur minimale du fond marin (m)
a   = 2000.0; % limite inferieur de l'onde en x (m)
b   = 5000.0; % limite superieur de l'onde en y (m)
Ly  = 2000.0; % longueur characteristique de l'onde en y (m)

% Parametres numeriques :
Nx = 64; 			 % nombre de node du maillage en x (extreme inclus)
Ny = 64; 			 % nombre de node du maillage en y (extreme inclus)

ComputeDt = false;   %ajout PERSO - option pour le calcul du dt
dt  = tfin/Nsteps;  %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 1.0;          % valeur de la condition CFL

type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite
mode_num_x = 2;			 % valeur propre en x
mode_num_y = 2;			 % valeur propre en y
bc_left   = 'harmonic';  % condition au bord gauche: dirichlet, neumann, harmonic
bc_right  = 'neumann';   % condition au bord droite: dirichlet, neumann, harmonic
bc_lower  = 'neumann';   % condition au bord inferieur: dirichlet, neumann, harmonic
bc_upper  = 'neumann';   % condition au bord superieur: dirichlet, neumann, harmonic
impulsion = false;		 % si vrais et harmonique une seule onde est emise
type_init = 'homogene';  % type initialisation: harmonic, default: homogene
F0 = 1.e0;			           % amplitude de l'harmonique initiale
A = 1.e0;			           % amplitude condition bord harmonique
omega = 5.e0;			       % frequence condition bord harmonique
write_mesh = true;		     % si vrais le maillage est ecrite dans output_mesh.out
write_f = true;		 	     % si vrais la solution est ecrite dans output_f.out
n_stride = 0;			     % nombre de stride

% Simulations
filename2  = "Nx_"+ num2str(Nx) + "Ny_" + num2str(Ny);
Nx_loc = Nx; Ny_loc = Ny; 
writeConfig;
disp('Exercice7_Kervyn_LeMeur configuration.in');   
system('Exercice7_Kervyn_LeMeur configuration.in'); 

%% Plot
ViewFormat; close all; facteurTemps = 4;
data = load(filename2+'_f.out');
mesh = load(filename2+'_mesh.out');
for i =1:length(data(:,1))/Ny
    surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
    title("t ="+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    zlim([-20 20]);
%     view(0,90);
    pause((1/facteurTemps)*dt);
end
