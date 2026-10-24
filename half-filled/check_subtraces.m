[U,UDIR] = QcmP.post.get_list('U')

for i=1:length(U)
   cd(UDIR(i))
   RDM12 = rdm_load('reduced_density_matrix_i1l1_i2l1.dat');
   RDM34 = rdm_load('reduced_density_matrix_i3l1_i4l1.dat');
   RDM13 = rdm_load('reduced_density_matrix_i1l1_i3l1.dat');
   RDM24 = rdm_load('reduced_density_matrix_i2l1_i4l1.dat');
   RDM14 = rdm_load('reduced_density_matrix_i1l1_i4l1.dat');
   RDM23 = rdm_load('reduced_density_matrix_i2l1_i3l1.dat');
   cd('..')

   % The unique index of each site (on the cartesian plane) is as follows:
   %
   % 3 4
   % 1 2

   fprintf('Checking subtraces for U=%.2f\n',U(i))
   if abs(norm(RDM12-RDM34))>=1e-4
       warning("We are breaking translational invariance along y (stripy phase?): RDM12 vs RDM34")
       fprintf('Difference norm for RDM12-RDM34: %.6e\n', abs(norm(RDM12-RDM34)))
   end
   if abs(norm(RDM13-RDM24))>=1e-4
       warning("We are breaking translational invariance along x (stripy phase?): RDM13 vs RDM24")
       fprintf('Difference norm for RDM13-RDM24: %.6e\n', abs(norm(RDM13-RDM24)))
   end
   if abs(norm(RDM12-RDM13))>=1e-4
       warning("We are breaking C4 symmetry (nematic phase?): RDM12 vs RDM13")
       fprintf('Difference norm for RDM12-RDM13: %.6e\n', abs(norm(RDM12-RDM13)))
   end
   if abs(norm(RDM34-RDM24))>=1e-4
       warning("We are breaking C4 symmetry (nematic phase?): RDM34 vs RDM24")
       fprintf('Difference norm for RDM34-RDM24: %.6e\n', abs(norm(RDM34-RDM24)))
   end
   if abs(norm(RDM14-RDM23))>=1e-4
       warning("We are breaking NNN symmetry (what?): RDM14 vs RDM23")
       fprintf('Difference norm for RDM14-RDM23: %.6e\n', abs(norm(RDM14-RDM23)))
   end

end

% contains

function rdm = rdm_load(fname)
   try
      rdm = QcmP.post.get_Hloc(fname);
   catch
      rdm = load(fname);
   end
end