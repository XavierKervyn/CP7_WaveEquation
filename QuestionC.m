%% Matlab - question c.
clear all; close all; clc;
Nsteps = 500;

% Parametres physiques:
tfin = 10;
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
pert_amplitude = 0.e0;	 pert_velocity  = 0.e0;	

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
omegaMN = pi*u*sqrt( (mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2 );%2*pi/T;

tfin = T;

% Simulations
filename2  = ['CA_Nx_',num2str(Nx),'Ny_',num2str(Ny)];
Nx_loc = Nx; Ny_loc = Ny; fname_loc = filename2;
writeConfig1;
disp('Exercice7_Kervyn_LeMeur configuration.in');   
system('Exercice7_Kervyn_LeMeur configuration.in'); 

% Plot simulation
ViewFormat;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);

figure('Name','plotsimulation')
for i =1:length(data(:,1))/Ny
    s = surfc(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]');
    colormap 'winter'
    zlim([-2 2]);
%     view(0,0);
    pause(1.e-5);
%     pause();
end
    SaveIMG("CSimulation");

figure('Name','plotAnalytique')
for i =1:length(data(:,1))/Ny
    surfc(mesh(1,:),mesh(2,:),F0*sin(pi*mode_num_x/(xR-xL)*mesh(1,:)).*sin(pi*mode_num_y/(yU-yL)*mesh(2,:)')*cos(omegaMN*data(Ny*i,1)));
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
    zlim([-2 2]);
%     view(0,0);
    pause(1.e-5);
%     pause();
end
    SaveIMG("CAnalytique");

%% soustraction des solutions
ViewFormat;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);
figure('Name','Soustraction')
for i =1:length(data(:,1))/Ny
    surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1)-F0*sin(pi*mode_num_x/(xR-xL)*mesh(1,:)).*sin(pi*mode_num_y/(yU-yL)*mesh(2,:)')*cos(omegaMN*data(Ny*i,1)));
    title("t ="+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
    colormap 'autumn'
%     zlim([-0.1 0.1]);
%     view(0,0);
    pause(1.e-9);
end
   SaveIMG("CdifferenceNumAna");

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

% Parametres numeriques :
% nombre de node du maillage en x,y (extreme inclus)
N = [16, 32, 64]; nsimul = 30;
Errors = zeros(length(N),nsimul);
BCFLs  = zeros(length(N),nsimul);
for n=1:length(N)
    Nx = N(n); 			 Ny = N(n);

    %période d'oscillation th.
    T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
    if(~ComputeDt)
        tfin = T;
        Nsteps = round(logspace(1.15,3,nsimul));
        dt = tfin./Nsteps; dt(1) = 1/u * sqrt(1/( ((Nx-1)/(xR-xL))^2 + ((Ny-1)/(yU-yL))^2 ));
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
    for k=1:nsimul
        data   = load(filename2(k)+'_f.out');
        mesh   = load(filename2(k)+'_mesh.out');
        %on récupère fnum(t=T), un tableau de taille (Nx)*Ny
        ErrF  = data(end-(Ny-1):end,2:end)-F0*sin(pi*mode_num_x/(xR-xL)*mesh(1,:)).*sin(pi*mode_num_y/(yU-yL)*mesh(2,:)').*cos(omegaMN*data(end-(Ny-1):end,1));
        for i = 1:Nx-1
            for j = 1:Ny-1
            Errors(n,k) = Errors(n,k) + ErrF(i,j)^2 + ErrF(i+1,j)^2 + ErrF(i,j+1)^2 + ErrF(i+1,j+1)^2;
            end
        end
    end
    Errors(n,:) = sqrt(0.25*((xR-xL)/(Nx-1))*((yU-yL)/(Ny-1))*Errors(n,:));
    BCFLs(n,:) = sqrt(u^2*dt.^2*(((Nx-1)/(xR-xL))^2 + ((Ny-1)/(yU-yL))^2));
end 

%%
ViewFormat;
figure('Name',"étude de convergence")
for n=1:length(N)
    plot(BCFLs(n,:),Errors(n,:), '+-','LineWidth',lw);
    hold on
end
    leg = legend(string(N),'Location','northwest');
    title(leg,"$N_{x,y}$")
    xlim([0,1]); ylim([0,0.015]);
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\beta_{CFL}$'); ylabel('$\varepsilon$ [a.u.]');
    SaveIMG("CconvergenceCFL2");
