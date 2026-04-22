function [ids,ens,X_list]  = energy_line(var,suffix,X_list)
%% Getting a list of energy values, from directories.
%
%     [ids,ens,X_list]  = QcmP.post.energy_line(var,suffix,X_list)
%
%  ids    : a cell of strings, the QcmPlab names for the pot-energy terms 
%  ens    : a cell of float-arrays, corresponding to the names above, forall X
%  var    : an optional char for the name of the line variable [default: 'U'] 
%  suffix : an optional charvec, handling inequivalent filename endings
%  X_list : an optional array of values for the line variable: where to search
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<1)
        var = 'U';
    end
    if(nargin<2)
        suffix = [];
    end
    if(nargin<3 || isempty(X_list))
        [X_list, ~] = QcmP.post.get_list(var); 
    else
        X_list = sort(X_list);
    end
    % Then we can proceed spanning all the U-values
    Nx = length(X_list);
    cellEn = cell(Nx,1);
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
        [ids, cellEn{i}] = QcmP.post.get_energies(suffix);
        cd('..');
    end
    % We need some proper reshaping
    Nids = length(ids);
    ens = cell(1,Nids);
    for jEn = 1:Nids
        ens{jEn} = zeros(Nx,1);
        for i = 1:Nx
           ens{jEn}(i) = cellEn{i}(jEn);
        end
        if(~isempty(suffix))
            filename = [ids{jEn},'_',suffix,'.txt'];
        else
            filename = [ids{jEn},'.txt'];
        end
        QcmP.post.writematrix(ens{jEn},filename,'Delimiter','tab');
    end
end

