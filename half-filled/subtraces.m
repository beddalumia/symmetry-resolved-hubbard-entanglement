[U,UDIR] = QcmP.post.get_list('U')

for i=1:length(U)
   cd(UDIR(i))
   ! sanitize_plaquette_subtraces_orbital_resolved_rdm | tee LOG.subtraces
   cd('..')
end