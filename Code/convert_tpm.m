function [ tpm ] = convert_tpm( tpm ,nnodes)
%   This function converts the 2^n x 2^n transition probability matrix into
%   the 2^n x n transition probability matrix in which the state of each
%   node is assumed to be independent of all other nodes.
    
    nstates = 2^nnodes;

    disp('before');
    disp(tpm);

    if size(tpm,2) == nstates   
        new_tpm = zeros(nstates,nnodes);
        for i = 1:nstates
            for j = 1:nnodes
                for k = 1:nstates
                    state = char(de2bi(k-1,nnodes)'+'0');
                    if strcmp(state(j),'1')
                        new_tpm(i,j) = new_tpm(i,j) + tpm(i,k);
                    end
                end
            end
        end    
        tpm = new_tpm;
    end
end

