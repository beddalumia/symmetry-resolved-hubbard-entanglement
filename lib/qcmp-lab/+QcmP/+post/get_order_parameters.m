function [names, order_parameters] = get_order_parameters()
%% Getting all information from order_parameters_[].dat
%
%       [names, order_parameters] = QcmP.post.get_order_parameters()
%
%  names: a cell of strings, the names of the order parameters (from [])
%  order_parameters: an array of floats, corresponding to the names above
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pattern = 'order_parameters_';
    files = dir('.');
    N = length(files);
    for i = 1:N
        temp_name = getfield(files,{i},'name');
        found = strfind(temp_name,pattern);
        if found
           full_name = temp_name;
        end
    end
    beheaded = erase(full_name,pattern); % Removes 'order_parameters_'
    detailed = erase(beheaded,'.dat');   % Removes '.dat'
    names = strsplit(detailed,'_');      % Reads the names separated by '_'
    order_parameters = load(full_name);
    if length(names) ~= length(order_parameters)
       error('Something went wrong reading order parameters!'); 
    end
end




