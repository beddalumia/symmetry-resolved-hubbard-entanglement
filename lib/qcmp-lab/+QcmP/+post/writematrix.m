function writematrix(variable,filename,varargin)
%% Wrapper of MATLAB builtin WRITEMATRIX, if not available (old version or
%  running on GNU Octave), it falls back to save() with --ascii option.
%
%           QcmP.post.writematrix(variable,filename,varargin)
%
%  EXAMPLE: QcmP.post.writematrix(Z,'zeta.txt','Delimiter','tab');
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    try % MATLAB >= R2019a
        writematrix(variable,filename,varargin{:});
    catch
        save(filename,'variable','-ascii','-double')
    end
end

