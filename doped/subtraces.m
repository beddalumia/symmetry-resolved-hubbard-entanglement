[mu,muDIR] = QcmP.post.get_list('xmu')

for i=1:length(mu)
   cd(muDIR(i))
   ! sanitize_plaquette_subtraces_orbital_resolved_rdm | tee LOG.subtraces
   cd('..')
end