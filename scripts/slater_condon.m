function smatrix = slater_condon(nmodes)
% slater_condon : Implementation of Slater-Condon rules for fermions
%                It pre-computes all the ❬istate| cdg_is c_js |jstate❭
%                matrix elements, storing them in a 4D array.
%                This can be used to build many-body representations of
%                one-body operators. NB: it assumes Nup,Ndw conservation.
%
%  >> smatrix = slater_condon(nmodes :: number of single-fermion modes)
%     smatrix :: 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]
%     
% ! This can be made way faster by implementing Slater-Condon
%   rules in an efficient machine-tuned way, as discussed in
%   https://arxiv.org/abs/1311.6244 (hal-01539072)
%
% © Gabriele Bellomia, 2025
    N = 4^nmodes;
    smatrix = zeros(2*nmodes,2*nmodes,N,N);
    for istate = 0:1:N-1 % nmode-orbital states
        for jstate = 0:1:N-1 % nmode-orbital states
            % ∑_s ∑_ij ❬istate| cdg_is c_js |jstate❭
            for ispin=1:2
                for imode=1:nmodes
                    for jmode=1:nmodes
                        % Apply cdg_is to ❬istate|
                        ibra = build_ket(istate,nmodes);
                        if ibra(imode+(ispin-1)*nmodes)==0
                            continue
                        end
                        ibra(imode+(ispin-1)*nmodes) = ibra(imode+(ispin-1)*nmodes)-1;
                        isign = (-1)^(sum(ibra(1:imode+(ispin-1)*nmodes)));
                        % Apply c_js to |jstate❭
                        jket = build_ket(jstate,nmodes);
                        if jket(jmode+(ispin-1)*nmodes)==0
                            continue
                        end
                        jket(jmode+(ispin-1)*nmodes) = jket(jmode+(ispin-1)*nmodes)-1;
                        jsign = (-1)^(sum(jket(1:jmode+(ispin-1)*nmodes)));
                        % Overlap ❬istate| cdg_is c_js |jstate❭
                        if isequal(ibra,jket)
                            smatrix(imode+(ispin-1)*nmodes,jmode+(ispin-1)*nmodes,istate+1,jstate+1) = isign * jsign;
                        end
                    end
                end
            end
        end
    end
end

function [vec,ket] = build_ket(state,Nlat)
  % BUILD_KET : Puts together a bit representation of a Slater determinant
  %             (and, if asked, a pretty string to print it to screen...)
  %  
  %  >> [vec,ket] = build_ket(state)
  %
  %  state :: integer representation of a basis state (bits are occupation numbers)
  %
  %  This depends entirely on the Fock basis conventions choosen in all
  %  ED-based solvers from QcmPlab (and many other codes, actually)
  %
  %  Copyright 2022 Gabriele Bellomia
  %
  Norb = 1;
  for ilat = 1:Nlat
      for ispin = 1:2
          shift = (ilat-1)*Norb + (ispin-1)*Norb*Nlat;
          index = shift+(1:Norb);
          vec(index) = bitget(state,index);
      end
  end
  if nargout>1
      kup = num2str(vec(1:Norb*Nlat));
      kdw = num2str(vec(Norb*Nlat+1:end));
      ket = ['| ',strrep(kup,'1','↑'),' 〉⊗ ',...
          '| ',strrep(kdw,'1','↓'), ' 〉'];
      ket = strrep(ket,'0','•');
  end
end

% FROM ED_SETUP:
% |imp_up>|bath_up> * |imp_dw>|bath_dw>        <- 2*Nlat*Norb bits
% |imp_sigma> = | (1…Norb)_1...(1…Norb)_Nlat > <--- Nlat*Norb bits
% lso indices are: io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
