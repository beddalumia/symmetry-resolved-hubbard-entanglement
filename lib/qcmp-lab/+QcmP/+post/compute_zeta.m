function [Z_list,X_list] = compute_zeta(filename,var,window,X_list)
%% Computing a list of Z-weight values, from self-energies on *real* axis.
%
%       [Z_list,X_list] = QcmP.post.compute_zeta(filename,window,X_list)
%
%  filename : an optional charvec, specifying which self-energy term to load
%  var      : an optional char for the name of the line variable [default: 'U'] 
%  window   : an optional window for the fit, abs(w)≤window [default 0.01]
%  Z_list   : a float-array, forall X, giving all the quasiparticle weight values
%  X_list   : an optional array of values for the line variable: where to search
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<2)
        var = 'U';
    end
    if ~exist('X_list','var') || isempty(X_list)
       [X_list, ~] = QcmP.post.get_list(var); 
    else
       X_list = sort(X_list);
    end
    if ~exist('window','var') || isempty(window)
       window = 0.01; 
    end
    % Then we can proceed spanning all the U-values
    Nx = length(X_list);
    Z_list = zeros(Nx,1);
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
        if(~exist('filename','var') || isempty(filename))
            allfiles = dir('impSigma*realw.ed');
            filename = allfiles(1).name;
        end
        if(~contains(filename,'Sigma'))
           warning('The file may not contain a self-energy.'); 
        end
        if(~contains(filename,'realw'))
           warning('This algorithm works for the real axis only.'); 
        end
        SIGMA = QcmP.plot.spectral_load(filename);
        Z_list(i) = fit_Zweight(SIGMA,window);
        cd('..');
    end
    suffix = erase(erase(filename,'impSigma_'),'realw.ed');
    if(~isempty(suffix))
        filename = ['zeta_',suffix,'realw.txt'];
    else
        filename = 'zeta_realw.txt';
    end
    QcmP.post.writematrix(Z_list,filename,'Delimiter','tab');
end

function Z = fit_Zweight(Sigma,w_th)
%% FIT_ZWEIGHT extracts the Landau quasiparticle weight from \Sigma(\omega)
%
%  Input:
%       Sigma : complex valued array, in the QcmP.plot.spectral_load format
%       w_th  : threshold for the fit, abs(w)≤w_th
%  Output:
%       Z     : 1 / ( 1 - d/dw[Re\Sigma]@w=0 )
%
%% BSD 3-Clause License
% 
%  Copyright (c) 2020, Gabriele Bellomia [adapted from MOTTlab]
%  All rights reserved.

    w = Sigma.zeta; % Extract frequencies
    xToFit = w(abs(w)<=w_th); % Define fit window (domain)
    yToFit = real(Sigma.real(abs(w)<=w_th)); % Define fit window (codomain)
    
    LinearModel = polyfit(xToFit,yToFit,1); % Fit to linear model
    
    %% Nozieres Theorem: Re{\Sigma(w)} = A + (1-1/Z)w + O(w^3)
    Z = 1/(1-LinearModel(1)); 

end