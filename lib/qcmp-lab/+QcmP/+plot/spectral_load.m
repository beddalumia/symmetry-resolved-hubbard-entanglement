function fstruct = spectral_load(filename)
%% SPECTRAL_LOAD: Loads a QcmPlab generated spectral tensor into named fields
%
%   >> fstruct = QcmP.plot.spectral_load(filename)
%
%  filename : filename of the complex spectral function to be plotted
%  fstruct  : data structure containing zeta(*), real and imag parts of file
%
%  (*) zeta is either real or matsubara frequency, depending on <filename>
%  ------------------------------------------------------------------------

    % Load complex function of real/imag frequency z, as a simple 3D array
    F = load(filename);

    % Assign domain, real and imag part according to QcmPlab convention
    fstruct.zeta = F(:,1);
    fstruct.real = F(:,3);
    fstruct.imag = F(:,2);

end
