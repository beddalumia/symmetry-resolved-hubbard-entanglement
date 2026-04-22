function [ids,obs,X_list] = observables_line(var,suffix,X_list)
%% Getting a list of observable values, from directories.
%
%       [ids,obs,X_list] = QcmP.post.observables_line(var,suffix,X_list)
%
%  ids    : a cell of strings, the QcmPlab names of the observables 
%  obs    : a cell of float-arrays, corresponding to the names above, forall X
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
    cellobs = cell(Nx,1);
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
        [ids, cellobs{i}] = QcmP.post.get_observables(suffix);
        cd('..');
    end
    % We need some proper reshaping
    Nobs = length(ids);
    obs = cell(1,Nobs);
    for jOBS = 1:Nobs
        obs{jOBS} = zeros(Nx,1);
        for i = 1:Nx
           obs{jOBS}(i) = cellobs{i}(jOBS);
        end
        if(~isempty(suffix))
            filename = [ids{jOBS},'_',suffix,'.txt'];
        else
            filename = [ids{jOBS},'.txt'];
        end
        QcmP.post.writematrix(obs{jOBS},filename,'Delimiter','tab');
    end
end

