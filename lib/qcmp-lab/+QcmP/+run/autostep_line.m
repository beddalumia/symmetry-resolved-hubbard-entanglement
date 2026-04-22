function autostep_line(EXE,doMPI,old,start,stop,var,varargin)
%% Runs a $(var)-line with feedback-controlled steps: 
%  if dmft does not converge the point is discarded and the step is reduced.
%
%   QcmP.run.autostep_line(EXE,doMPI,old,start,stop,var,varargin)
%
%   EXE                 : Executable driver
%   doMPI               : Flag to activate OpenMPI
%   old                 : Restart point [NaN or empty -> no restart]
%   start,stop          : Input values for $(var) [start<x<stop or start>x>stop]
%   var                 : Main variable of the line [default: 'U']
%   varargin            : Set of fixed control parameters ['name',value]

%% Print docstring if no input is provided
if nargin < 1
   help QcmP.run.autostep_line 
   return
end

xlist = fopen(sprintf('%s_list.txt',var),'a');
xconv = fopen(sprintf('%s_conv.txt',var),'a');

%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step = [.5,.25,.1,.05,.01,.005,.001];  % Let's keep them hard-coded... for now.
Nstep = length(step);

% Autodetermine span direction
if start > stop
   step = -step;
end

nonconvFLG = false;                % Convergence-fail *flag*
nonconvCNT = 0;                    % Convergence-fail *counter*
nonconvMAX = Nstep-1;              % Maximum #{times} we accept DMFT to fail

x = start; 
while abs(x-step-stop) > min(abs(step)/2)
%                  ≥ 0 would sometimes give precision problems

    unconverged = QcmP.run.single_point(EXE,doMPI,var,x,old,varargin{:});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
    if (unconverged) 
        nonconvFLG = true;
        nonconvCNT = nonconvCNT + 1;
        errorfile = [sprintf('%s=%f',var,x),'/ERROR.README'];
        movefile(errorfile,sprintf('ERROR_%s=%f',var,x));
        fprintf(xlist,'%f\n', NaN);             % Write on $(var)_list
    else
        fprintf(xlist,'%f\n', x);	            % Write on $(var)-list
        fprintf(xconv,'%f\n', x);	            % Write on $(var)-conv
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nonconvCNT > nonconvMAX
        error('Not converged: phase-span stops now!');         
    end 

    if nonconvFLG == true
        x = old; 			        % if nonconverged we don't want to update 
        nonconvFLG = false;	        % > but we want to reset the flag(!)
    else
        old = x; 			        % else we update old and proceed to...         
    end

    x = x + step(nonconvCNT+1);     % ...update the line variable 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(xlist); fclose(xconv);

end





