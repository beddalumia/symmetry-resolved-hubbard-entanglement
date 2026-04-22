function [ids,ordpms,X_list] = order_parameter_line(var,X_list)
%% Getting a list of order parameter values, from directories.
%
%       [ids,ordpms,X_list] = QcmP.post.order_parameter_line(var,X_list)
%
%  ids    : a cell of strings, the names of the order parameters 
%  ordpms : a cell of arrays, corresponding to the names above, for all X
%  var    : an optional char for the name of the line variable [default: 'U'] 
%  X_list : an optional array of values for the line variable: where to search
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<1)
        var = 'U';
    end
    if(nargin<2 || isempty(X_list))
        [X_list, ~] = QcmP.post.get_list('U'); 
    else
        X_list = sort(X_list);
    end
    % Then we can proceed spanning all the U-values
    Nx = length(X_list);
    cellordpms = cell(Nx,1);
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
        [ids, cellordpms{i}] = QcmP.post.get_order_parameters();
        cd('..');
    end
    % We need some proper reshaping
    Nordpms = length(ids);
    ordpms = cell(1,Nordpms);
    for jORDPMS = 1:Nordpms
        ordpms{jORDPMS} = zeros(Nx,1);
        for i = 1:Nx
           ordpms{jORDPMS}(i) = cellordpms{i}(jORDPMS);
        end
        filename = [ids{jORDPMS},'.txt'];
        QcmP.post.writematrix(ordpms{jORDPMS},filename,'Delimiter','tab');
    end
end


