function [kins,X_list]  = kinetic_line(var,X_list)
%% Getting a list of kinetic energy values, from directories.
%
%     [kins,X_list] = QcmP.post.kinetic_line(var,X_list)
%
%  var    : an optional char for the name of the line variable [default: 'U'] 
%  kins   : an array of values for Kinetic energies, forall X
%  X_list : an optional array of values for the line variable: where to search
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<1)
        var = 'U';
    end
    if ~exist('X_list','var') || isempty(X_list)
       [X_list, ~] = QcmP.post.get_list(var); 
    else
       X_list = sort(X_list);
    end
    % Then we can proceed spanning all the U-values
    Nx = length(X_list);
    kins = zeros(Nx,1);
    for i = 1:length(X_list)
        x = X_list(i);
        DIR= sprintf('%s=%f',var,x);
        if ~isfolder(DIR)
            errstr = sprintf('%s_list appears to be inconsistent: ',var);
            errstr = [errstr,DIR]; %#ok
            errstr = [errstr,' folder has not been found.']; %#ok
            error(errstr);
        end
        cd(DIR); 
        if not(isfile('dmft_kinetic_energy.dat'))
           kins(i) = NaN;
           cd('..');
           continue;
        end
        % The dmft_kinetic_energy.dat file has a weird structure...
        system('head -1 dmft_kinetic_energy.dat > tmp_energia_cinetica.txt');
        tempVec = load('tmp_energia_cinetica.txt'); % this should be safe...
        kins(i) = tempVec(1); % For sure the 1st value is the total K.E.
        system('rm tmp_energia_cinetica.txt');
        cd('..');
    end
    filename = 'kinetic_energy.txt';
    QcmP.post.writematrix(kins,filename,'Delimiter','tab');
end




