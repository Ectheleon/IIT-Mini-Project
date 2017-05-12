function [ dat ] = genSeries( Trans_Matrix, n )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    state = 0;
    
    nsettle = 500;
    
    dat = zeros([n 3]);
    
    tmp = cumsum(Trans_Matrix,2);
    
    for i=1:nsettle+n
        p = rand(1);
        ind = 1;
        while tmp(state+1, ind)<p
            ind = ind+1;
        end
        
        state = ind-1;
        
        if i>nsettle
            for bit=1:3
                dat(i-nsettle,bit) = bitget(state,bit);
            end
        end
        
    end


end

