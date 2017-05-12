k = [0.1 0.1 0.1];
p = [0.1 0.1];


%
model = 1;
for i =1:9
    
    
    for j = 1:9
        p(1)= 0.1*j;
        p(2) = k(1)^2 / (1.0*p(1));
        if i~=j
            gen_mat(model, p,k,1);
        end 
    end
    k = k + 0.1;

end

%%

num = 0.01;
k = [num num num];
p = [num num];

show = gen_mat(3, p, k,1);