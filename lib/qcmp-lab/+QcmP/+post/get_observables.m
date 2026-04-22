function [names, observables] = get_observables(suffix)
%% Getting all information from observables_last.ed and observables_info.ed
%
%       [names, observables] = QcmP.post.get_observables(suffix)
%
%  names: a cell of strings, the QcmPlab names of the observables
%  observables: an array of float values, corresponding to the names above
%  suffix: an optional charvec string, handling inequivalent filename endings
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(~exist('suffix','var'))
        suffix = [];
    end
    if(~isempty(suffix))
        filename = ['observables_last_',suffix,'.ed'];
    else
        filename = 'observables_last.ed';
    end
    observables = load(filename);
    try % MATLAB >= R2019a
        names = readcell('observables_info.ed','FileType','fixedwidth');
        names(strcmp(names,'#'))=[];
        for i = 1:length(names)
            tempstr = names{i};                 % Temporary string variable
            head = sscanf(tempstr,'%d');        % Extracts the initial integer
            head = int2str(head);               % Int to Str conversion
            names{i} = erase(tempstr,head);     % Proper beheading ;D
        end
    catch 
        names = cell(length(observables),1);
        for i=1:length(observables)
            names{i} = sprintf('observable#%d',i);
        end
        warning('Could not read observables_info.ed, you might want to look at it yourself.')
    end
end





