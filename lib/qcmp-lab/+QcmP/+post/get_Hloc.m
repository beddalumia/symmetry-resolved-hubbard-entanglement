function H = get_Hloc(filename)
%% GET_HLOC reads a complex matrix in the "TB_write_Hloc" format
%
%   >> H = get_Hloc(<filename>)
%
%   filename: optional string for custom filenames (default="Hloc.txt")
    if nargin < 1
        filename = "Hloc.txt";
    end 
    F = load(filename);
    N = min(size(F));
    H = F(1:N,:);
    H = H + 1i*F(N+1:end,:);
end

