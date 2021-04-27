%% Matlab - question c.
clear all; close all; clc;
Nsteps = 1000;

% Parametres physiques:
tfin = 20.0;    % temps physique de la simulation (s)
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
pert_amplitude = 3.e0;	 pert_velocity  = 0.e0;	

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
% nombre de node du maillage en x,y (extreme inclus)
Nx = 64; 			 Ny = 64;

ComputeDt = false;   %ajout PERSO - option pour le calcul du dt
dt  = tfin/Nsteps;  %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 1.0;          % valeur de la condition CFL

type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite

% valeur propre en x,y
mode_num_x = 3;			 mode_num_y = 5;			 

%conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'dirichlet';  bc_right  = 'dirichlet';   
bc_lower  = 'dirichlet';  bc_upper  = 'dirichlet';   

impulsion = false;          % si vrais et harmonique une seule onde est emise
type_init = 'harmonic';     % type initialisation: harmonic, default: homogene
F0 = 1.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 10.e0;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrais le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrais la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

%période d'oscillation th.
T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
nbperiodes = 5;
if(ComputeDt)
    tfin = nbperiodes*T;
    dt = tfin/Nsteps;
end
omegaMN = 2*pi/T;

% Simulations
filename2  = ['CA_Nx_',num2str(Nx),'Ny_',num2str(Ny)];
Nx_loc = Nx; Ny_loc = Ny; fname_loc = filename2;
writeConfig1;
disp('Exercice7_Kervyn_LeMeur configuration.in');   
system('Exercice7_Kervyn_LeMeur configuration.in'); 

% Plot simulation
ViewFormat; facteurTemps = 8;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);

figure('Name','plotsimulation')
for i =1:length(data(:,1))/Ny
    surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
    zlim([-5 5]);
%     view(0,0);
    pause((1/facteurTemps)*dt);
%     pause();
end

figure('Name','plotAnalytique')
[X, Y]  = meshgrid(mesh(1,:),mesh(2,:));
for i =1:length(data(:,1))/Ny
    surf(X,Y,F0*sin(pi*mode_num_x/(xR-xL)*X).*sin(pi*mode_num_y/(yU-yL)*Y)*cos(omegaMN*data(Ny*i,1)));
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
    zlim([-5 5]);
%     view(0,90);
    pause((1/facteurTemps)*dt);
end

% soustraction des solutions
ViewFormat; facteurTemps = 8;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);
figure('Name','Soustraction')
[X, Y]  = meshgrid(mesh(1,:),mesh(2,:));
for i =1:length(data(:,1))/Ny
    surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1)-F0*sin(pi*mode_num_x/(xR-xL)*X).*sin(pi*mode_num_y/(yU-yL)*Y)*cos(omegaMN*data(Ny*i,1)));
    title("t ="+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
%     zlim([-0.1 0.1]);
%     view(0,0);
    pause((1/facteurTemps)*dt);
%     pause();
end

%% Matlab - question c: étude de convergence
clear all; close all; clc;

% Parametres physiques:
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
pert_amplitude = 0.e0;	 pert_velocity  = 0.e0;	

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
Nx = 16; 			 Ny = 16;

ComputeDt = false;   %ajout PERSO - option pour le calcul du dt
CFL = 1.0;          % valeur de la condition CFL

type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite

% valeur propre en x,y
mode_num_x = 2;			 mode_num_y = 2;			 

%conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'dirichlet';  bc_right  = 'dirichlet';   
bc_lower  = 'dirichlet';  bc_upper  = 'dirichlet';   

impulsion = false;          % si vrais et harmonique une seule onde est emise
type_init = 'harmonic';     % type initialisation: harmonic, default: homogene
F0 = 1.e0;			        % amplitude de l'harmonique initiale
A  = 1.e0;			        % amplitude condition bord harmonique
omega = 10.e0;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrai le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrai la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

%période d'oscillation th.
T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
if(~ComputeDt)
    tfin = T;
    nsimul = 30; Nsteps = round(logspace(2,3,nsimul));
    dt = tfin./Nsteps;
end
omegaMN = 2*pi/T;

% Simulations
filename2  = "Nsteps_"+ string(Nsteps);

for i=1:nsimul
    N_loc  = Nsteps(i);
    dt_loc = dt(i);
    writeConfig2;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 
end

% calcul de l'erreur
NperErr = zeros(1,nsimul);
for k=1:nsimul
    data   = load(filename2(k)+'_f.out');
    mesh   = load(filename2(k)+'_mesh.out');
    [X, Y] = meshgrid(mesh(1,:),mesh(2,:));
    %on récupère fnum(t=T), un tableau de taille (Nx)*Ny
    ErrF  = data(end-(Ny-1):end,2:end)-F0*sin(pi*mode_num_x/(xR-xL)*X).*sin(pi*mode_num_y/(yU-yL)*Y).*cos(omegaMN*data(end-(Ny-1):end,1));
    for i = 1:Nx-1
      for j = 1:Ny-1
        NperErr(k) = NperErr(k) + ErrF(i,j)^2 + ErrF(i+1,j)^2 + ErrF(i,j+1)^2 + ErrF(i+1,j+1)^2;
      end
    end
end
NperErr = sqrt(0.25*((xR-xL)/(Nx-1))*((yU-yL)/(Ny-1))*NperErr); 

%
ViewFormat;
figure('Name',"étude de convergence")
    plot(1./(Nsteps.^2),NperErr, '+-','LineWidth',lw);
    grid minor; set(gca,'fontsize',fs);
    xlabel('$1/(N_{steps})^2$'); ylabel('$\varepsilon$ [a.u.]');     
