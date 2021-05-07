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
ViewFormat; levels = 15;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);

% figure('Name','plotsimulation')
[X, Y]  = meshgrid(mesh(1,:),mesh(2,:));
axes1 = axes('Parent',figure);
for i =1:length(data(:,1))/Ny
    s1 = contourf(X,Y,data(Ny*(i-1)+1:Ny*i,2:Nx+1),levels);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); 
    title("$t =$ "+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    set(gca,'fontsize',fs)
    caxis([-(A+0.2) A+0.2]);
        h = colorbar(axes1);
        h.Label.String = '$f$ [a.u.]';
        h.Label.Interpreter = 'Latex';
        h.Label.FontSize = fs;
end
    SaveIMG("CSimulation");

% figure('Name','plotAnalytique')
[X, Y]  = meshgrid(mesh(1,:),mesh(2,:),15);
axes1 = axes('Parent',figure);
for i =1:length(data(:,1))/Ny
    s1 = contourf(X,Y,F0*sin(pi*mode_num_x/(xR-xL)*mesh(1,:)).*sin(pi*mode_num_y/(yU-yL)*mesh(2,:)')*cos(omegaMN*data(Ny*i,1)),levels);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); 
    title("$t =$ "+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    set(gca,'fontsize',fs)
    caxis([-(A+0.2) A+0.2]);
        h = colorbar(axes1);
        h.Label.String = '$f$ [a.u.]';
        h.Label.Interpreter = 'Latex';
        h.Label.FontSize = fs;
end
    SaveIMG("CAnalytique");

%% soustraction des solutions
ViewFormat; levels = 15;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);
[X, Y]  = meshgrid(mesh(1,:),mesh(2,:));
axes1 = axes('Parent',figure);
for i =1:length(data(:,1))/Ny
    s1 = contourf(X,Y,data(Ny*(i-1)+1:Ny*i,2:Nx+1)-F0*sin(pi*mode_num_x/(xR-xL)*mesh(1,:)).*sin(pi*mode_num_y/(yU-yL)*mesh(2,:)')*cos(omegaMN*data(Ny*i,1)),levels);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); 
    title("$t =$ "+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    set(gca,'fontsize',fs)
        h = colorbar(axes1);
        h.Label.String = '$f_{num} - f_{ana}$ [a.u.]';
        h.Label.Interpreter = 'Latex';
        h.Label.FontSize = fs;
    colormap 'autumn'
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

for i=1:length(BCFLs(3,:))
    if(BCFLs(2,i) > 1)
       Errors(2,i) = 0; 
    end
    if(BCFLs(3,i) > 1)
       Errors(3,i) = 0; 
    end
end
%%
ViewFormat;
figure('Name',"étude de convergence")
for n=1:length(N)
    plot(BCFLs(n,:),Errors(n,:), '+-','LineWidth',lw);
    hold on
end
    leg = legend(string(N),'Location','northwest');
    title(leg,"$n_{x,y}$")
    xlim([0,1]); ylim([0,0.015]);
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\beta_{CFL}$'); ylabel('$\varepsilon$ [(a.u.)$\cdot$ m]');
    SaveIMG("CconvergenceCFL2");
%%
   BetaZoom1 = BCFLs(1,8:17); ErrZoom1 = Errors(1,8:17);
   BetaZoom2 = BCFLs(2,13:26); ErrZoom2 = Errors(2,13:26);
   BetaZoom3 = BCFLs(3,18:end); ErrZoom3 = Errors(3,18:end);
   
figure('Name',"convergence")
    plot(BetaZoom1,ErrZoom1, '+','LineWidth',lw);
    hold on
    plot(BetaZoom2,ErrZoom2, '+','LineWidth',lw);
    plot(BetaZoom3,ErrZoom3, '+','LineWidth',lw);
    
    P1 = polyfit(BetaZoom1,ErrZoom1,1); 
    z1  = polyval(P1, BetaZoom1);
    plot(BetaZoom1,z1,'--','Linewidth',lw,'HandleVisibility','off');
    
    P2 = polyfit(BetaZoom2,ErrZoom2,1); 
    z2  = polyval(P2, BetaZoom2);
    plot(BetaZoom2,z2,'--','Linewidth',lw,'HandleVisibility','off');

    P3 = polyfit(BetaZoom3,ErrZoom3,1); 
    z3  = polyval(P3, BetaZoom3);
    plot(BetaZoom3,z3,'--','Linewidth',lw,'HandleVisibility','off');
    
    leg = legend(string(N),'Location','northwest');
    title(leg,"$n_{x,y}$")
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\beta_{CFL}$'); ylabel('$\varepsilon$ [(a.u.)$\cdot$ m]');
SaveIMG("CconvergenceCFL2_lines");
   