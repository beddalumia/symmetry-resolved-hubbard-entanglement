function [states, probs] = get_pure_states(suffix,fast)
%% Getting all information about pure states (cross-checking with the RDM)
%
%       [states, probs] = QcmP.post.get_pure_states(suffix[,check])
%
%  states: a cell of pure-state vectors [states{i} is a complex array]
%  probs: an array of pure-state probabilities [probs(i) = 𝓟(states{i})]
%  suffix: a *required* charvec / string, handling inequivalent filenames
%  fast: an optional boolean to avoid RDM diagonalization [default:F]
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<1)
        help QcmP.post.get_pure_states
        return
    end
    if(nargin<2)
        fast = false; 
    end
    if(fast)
        if(~isempty(suffix))
            filename_p = ['probabilities_',suffix,'.dat'];
            filename_s = ['pure-states_',suffix,'.dat'];
        else
            filename_p = 'probabilities.dat';
            filename_s = 'pure-states.dat';
        end
        probs = load(filename_p);
        psmat = load(filename_s);
        Npure = length(probs);
    else
        filename_rdm = ['reduced_density_matrix_',suffix,'.dat'];
        try
            RDM = QcmP.post.get_Hloc(filename_rdm);
        catch
            RDM = load(filename_rdm); % imag part may be omitted if zero
        end
        Npure = min(size(RDM));
        %% ACHTUNG: our convention [probs(i) = 𝓟(states{i})] implies that
        %%          the matrix of states does not follow the convention of
        %%          linear algebra, i.e. usually the eigenvector matrix is
        %%          defined such that A=V*D*V', but for us A=V'*D*V [!🚨!]
        [psmat,probs] = eig(RDM,'vector'); psmat = psmat';
    end
    %% Reshape to cell
    states = cell(Npure,1);
    for i = 1:Npure
       states{i} = psmat(i,:);
    end
end


