function DessineSurface(nsimul) %% simule sur le Meshgrid
N1_grid = linspace(100,1000,nsimul); 
N2_grid = linspace(100,1000,nsimul); 
[X, Y]  = meshgrid(N1_grid,N2_grid);

N1_ = N1_grid; N2_ = N2_grid;
ViewFormat;

trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0;
p_          = 1.e0;
propMesh_   = false;

filename2 = strings;
for i = 1:nsimul
    for j=1:nsimul
        filename2(i,j) = "N1_"+ num2str(N1_(i))+"N2_"+num2str(N2_(j)) ;
    end
end

% Simulations 2 (non trivial)
for i = 1:nsimul
    for j=1:nsimul
        N1_loc = N1_(i);
        N2_loc = N2_(j);
        writeConfig;
        disp('Exercice6 configuration.in');   
        system('Exercice6 configuration.in');
    end
end

figure('Name','surface')
phiAna = 183.4148;
phirb  = zeros(nsimul,nsimul); %valeur de phi en r = b
for i=1:nsimul
    for j=1:nsimul
        data = load(filename2(i,j)+'_phi.out');
        r    = data(:,1);
        [val,indice] = min(abs(r - b_));
        phirb(i,j)     = abs(data(indice,2)-phiAna)/phiAna * 100;
    end
end
    s1 = surf(X,Y,phirb);
    xlabel('$N_2$');
    ylabel('$N_1$');
    zlabel('$\chi_b$ [\%]');
    colormap autumn
    view(-130,35);
    s1.FaceColor = 'interp';
    set(gca,'fontsize',fs)
SaveIMG("SurfaceN1N2");
end