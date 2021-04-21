%% Matlab: question e 
clear; close all; clc;
nsimul=1;
% Parametres physiques:
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
pert_amplitude = 3.e1;	 pert_velocity  = linspace(1,11.1,nsimul);	

% Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

% Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 2.e0;      % profondeur maximale du fond marin (m)
h1  = 1.e0;      % profondeur minimale du fond marin (m)
a   = 2000.0;    % limite inferieur de l'onde en x (m)
b   = 5000.0;    % limite superieur de l'onde en y (m)
Ly  = 2000.0;    % longueur characteristique de l'onde en y (m)

% Parametres numeriques :
% nombre de node du maillage en x,y (extreme inclus)
Nx = 64; 			 Ny = 64;

ComputeDt = false;                   %ajout PERSO - option pour le calcul du dt
CFL = 1;         % valeur de la condition CFL


type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite

% valeur propre en x,y
mode_num_x = 3;			 mode_num_y = 5;			 

%conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'dirichlet';  bc_right  = 'dirichlet';   
bc_lower  = 'dirichlet';  bc_upper  = 'dirichlet';   

impulsion = false;          % si vrais et harmonique une seule onde est emise
type_init = 'homogene';     % type initialisation: harmonic, default: homogene
F0 = 1.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 10.e0;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrai le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrai la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

%période d'oscillation th.
T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
nbperiodes = 20; tfin = nbperiodes*T;
omegaMN = 2*pi/T;
Nsteps = 2000; dt = tfin/Nsteps;

% Simulations
filename2  = "PVelocity_"+ string(pert_velocity);

for i=1:nsimul
    PV_loc = pert_velocity(i);
    writeConfig4;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 
end

%
ViewFormat; facteurTemps = 8;
for k=1:nsimul
    data = load(filename2(k)+'_f.out');
    mesh = load(filename2(k)+'_mesh.out');
    figure('Name',"plotsimulation PVelocity="+num2str(pert_velocity(k)))
    for i =1:length(data(:,1))/Ny
        surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
        title("$t =$"+num2str(data(Ny*i,1))+"s, $PV =$"+num2str(CFL(k))+"m/s",'Interpreter','latex');
        grid minor; set(gca,'fontsize',fs);
        xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
        zlim([-0.5 0.5]);
    %     view(0,0);
        pause((1/facteurTemps)*dt);
    %     pause();
    end
end