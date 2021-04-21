%va écrire le configuration.in
filename = 'configuration.in';
fid = fopen(filename,'wt'); 

%Paramètres physiques:
fprintf(fid,['tfin=', num2str(tfin,8),'\n']);
fprintf(fid,['xL=', num2str(xL),'\n']);
fprintf(fid,['xR=', num2str(xR),'\n']); 
fprintf(fid,['yL=', num2str(yL),'\n']); 
fprintf(fid,['yU=', num2str(yU),'\n']); 

fprintf(fid,['pert_amplitude=', num2str(pert_amplitude),'\n']);
fprintf(fid,['pert_velocity=', num2str(pert_velocity),'\n']);

%Paramètres physiques milieu uniforme
fprintf(fid,['u=', num2str(u),'\n']);

%Paramètres physiques onde Belharra
fprintf(fid,['g=', num2str(g),'\n']);
fprintf(fid,['h0=', num2str(h0),'\n']);
fprintf(fid,['h1=', num2str(h1),'\n']);
fprintf(fid,['a=', num2str(a),'\n']);
fprintf(fid,['b=', num2str(b),'\n']);
fprintf(fid,['Ly=', num2str(Ly),'\n']);

%Paramètres numériques:
fprintf(fid,['Nx=', num2str(Nx),'\n']);
fprintf(fid,['Ny=', num2str(Ny),'\n']);

fprintf(fid,['ComputeDt=', num2str(ComputeDt),'\n']);
fprintf(fid,['dt=', num2str(dt_loc,8),'\n']);
fprintf(fid,['CFL=', num2str(CFL_loc,8),'\n']);

fprintf(fid,['type_u2=', type_u2,'\n']);
fprintf(fid,['ecrire_f=', num2str(ecrire_f),'\n']);
fprintf(fid,['mode_num_x=', num2str(mode_num_x),'\n']);
fprintf(fid,['mode_num_y=', num2str(mode_num_y),'\n']);

fprintf(fid,['bc_left=', bc_left,'\n']);
fprintf(fid,['bc_right=',bc_right,'\n']);
fprintf(fid,['bc_lower=',bc_lower,'\n']);
fprintf(fid,['bc_upper=',bc_upper,'\n']);

fprintf(fid,['impulsion=', num2str(impulsion),'\n']);
fprintf(fid,['type_init=', type_init,'\n']);
fprintf(fid,['F0=', num2str(F0),'\n']);
fprintf(fid,['A=', num2str(A),'\n']);
fprintf(fid,['omega=', num2str(omega),'\n']);
fprintf(fid,['write_mesh=', num2str(write_mesh),'\n']);
fprintf(fid,['write_f=', num2str(write_f),'\n']);
fprintf(fid,['n_stride=', num2str(n_stride),'\n']);

%écriture dans des fichiers:
fprintf(fid,['output_mesh= CFL_',num2str(CFL_loc),'_mesh.out\n']); 
fprintf(fid,['output_file= CFL_',num2str(CFL_loc),'_f.out\n']); 
fprintf(fid,['output_energy= CFL_',num2str(CFL_loc),'_E.out\n']); 
fprintf(fid,['output_velocity= CFL_',num2str(CFL_loc),'_u.out\n']);

fclose(fid);