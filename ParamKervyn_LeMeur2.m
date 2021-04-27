%% PARTIE 7.3
% Vague avec profondeur de l'océan variable.
% question a.
clear all; close all; clc;
Nsteps = 1000;

% Parametres physiques:
tfin = 100.0;    % temps physique de la simulation (s)
xL   = 0.0 ; 		     xR   = 5.e3 ; 		     
yL   = 0.0 ;		     yU   = 2.e3;		%dimensions (m)     
pert_amplitude = 0.e0;	 pert_velocity  = 0.e0;	

% Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

% Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 1.e3;      % profondeur maximale du fond marin (m)
h1  = 2.e1;      % profondeur minimale du fond marin (m)
a   = 2000.0; % limite inferieur de l'onde en x (m)
b   = 5000.0; % limite superieur de l'onde en y (m)
Ly  = 2000.0; % longueur characteristique de l'onde en y (m)

% Parametres numeriques :
% nombre de node du maillage en x,y (extreme inclus)
Nx = 64; 			 Ny = 64;

ComputeDt = false;  %ajout PERSO - option pour le calcul du dt
dt  = tfin/Nsteps; %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 0.95;        % valeur de la condition CFL

type_u2  = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;         % if true f est ecrite

% valeur propre en x,y
mode_num_x = 0;			 mode_num_y = 0;			 

%conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'harmonic'; bc_right  = 'neumann';   
bc_lower  = 'neumann';  bc_upper  = 'neumann';   

T = 15; %2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));

impulsion = false;          % si vrais et harmonique une seule onde est emise
type_init = 'homogene';     % type initialisation: harmonic, default: homogene
F0 = 0.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 2*pi/T;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrais le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrais la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

% Simulations
filename2  = "Nx_"+ num2str(Nx) + "Ny_" + num2str(Ny);
    Nx_loc = Nx; Ny_loc = Ny; 
    writeConfig;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 

% Plot simulation
ViewFormat;
data = load(filename2+'_f.out');
mesh = load(filename2+'_mesh.out');
[X,Y]= meshgrid(mesh(1,:),mesh(2,:));

% profil de la profondeur
hK1 = h0*ones(length(X),length(Y));
for j = 1:length(Y)
    for i=1:length(X)
        x = mesh(1,i); %disp(x);
        if(x>a && x<b) 
            hK1(j,i) = h0+(h1-h0)*sin(pi*(x-a)/(b-a));
        end
    end
end

figure('Name','plotProfondeur Cas 1')
for i =1:length(data(:,1))/Ny
    surf(X,Y,hK1);
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
%     zlim([-5 5]);
%     view(0,0);
end

facteurTemps = 8;
figure('Name','plotsimulation')
for i =1:length(data(:,1))/Ny
    surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [m]')
    zlim([-5, 5]);   % view(0,0);
    pause((1/facteurTemps)*dt);
%     pause();
end