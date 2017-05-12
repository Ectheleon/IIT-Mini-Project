function [ trans_mat] = shuffle( k )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    trans_mat = ones(8);


    for i=0:7
        for j = 0:7
            for bit = 1:3
                if bitand(i,j)+bitxor(i,j) ==i
                    if bitget(i,bit)
                        if ~bitget(j,bit)
                            trans_mat(i+1, j+1) = trans_mat(i+1, j+1)*k(bit);
                        else
                            trans_mat(i+1, j+1) = trans_mat(i+1, j+1)*(1-k(bit));
                        end
                    end
                else
                    trans_mat(i+1, j+1) = 0;
                end
            end
        end
    end
end

