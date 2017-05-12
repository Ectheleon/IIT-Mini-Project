function [ tpm ] = gen_mat( model, p, k ,print)
        
    tpm = zeros(8);
    temp = zeros([8 4]);
    
    pmat = shuffle(k);

    Prob = [(1-p(1))*(1-p(2)) p(1)*(1-p(2)) p(2)*(1-p(1)) p(1)*p(2)];
    
    switch model
%%
% The first model deals with the case where workers are only assigned if
% there is a perfect fit.
        case 1 
            
            temp(:,1) = 0:7;
            
            for i=0:7
                for P = 1:3
                    temp(i+1, P+1) = temp(i+1,1);

                    if ~bitget(i,3)
                        if P==1 
                            temp(i+1,P+1) = temp(i+1,P+1)+4;
                        end
                    end

                    if ~bitget(i,2)
                        if P==2 
                            temp(i+1,P+1) = temp(i+1,P+1) + 2;
                        end
                    end

                    if ~bitget(i,1)
                        if P==3
                            temp(i+1,P+1) = temp(i+1,P+1) + 1;
                        end
                    end
                end
            end

            for i=1:8
                for j = 1:8
                    for Pind=1:4
                        ind = int8(temp(j,Pind))+1;
                        tpm(i,ind)=tpm(i,ind)+pmat(i,j)*Prob(Pind);
                    end
                end
            end
 
 %%
 % The second model deals with the case where workers are assigned to all
 % projects that can be completed without overflow. Thus two specialists
 % may be assigned to work together on a project, but a multi-skilled
 % person wont have to do a specialist task.
        case 2

            temp(:,1) = 0:7;

            for i=0:7
                for P = 1:3
                    temp(i+1, P+1) = temp(i+1,1);

                    if ~bitget(i,3)
                        if P==1 || (P==3 && bitget(i,1) && ~bitget(i,2))
                            temp(i+1,P+1) = temp(i+1,P+1)+4;
                        end
                    end

                    if ~bitget(i,2)
                        if P==2 || (P==3 && bitget(i,1) &&~bitget(i,3))
                            temp(i+1,P+1) = temp(i+1,P+1) + 2;
                        end
                    end

                    if ~bitget(i,1)
                        if P==3
                            temp(i+1,P+1) = temp(i+1,P+1) + 1;
                        end
                    end
                end
            end

            for i=1:8
                for j = 1:8
                    for Pind=1:4
                        ind = int8(temp(j,Pind))+1;
                        tpm(i,ind)=tpm(i,ind)+pmat(i,j)*Prob(Pind);
                    end
                end
            end
 %%
 % This model deals with the case where workers are always assigned to
 % projects if there is a way to complete the project. Waste is not
 % considered a problem if the proeject can be completed.     
        case 3
            
            temp(:,1) = 0:7;
            
            for i=0:7
                for P = 1:3
                    temp(i+1, P+1) = temp(i+1,1);

                    if ~bitget(i,3)
                        if P==1 || (P==3 && bitget(i,1) && ~bitget(i,2))
                            temp(i+1,P+1) = temp(i+1,P+1)+4;
                        end
                    end

                    if ~bitget(i,2)
                        if P==2 || (P==3 && bitget(i,1) &&~bitget(i,3))
                            temp(i+1,P+1) = temp(i+1,P+1) + 2;
                        end
                    end

                    if ~bitget(i,1)
                        if P==3 || (bitget(i,2) && bitget(i,3))
                            temp(i+1,P+1) = temp(i+1,P+1) + 1;
                        end
                    end
                end
            end


            for i=1:8
                for j = 1:8
                    for Pind=1:4
                        ind = int8(temp(j,Pind))+1;
                        tpm(i,ind)=tpm(i,ind)+pmat(i,j)*Prob(Pind);
                    end
                end
            end           
%%
% In all other cases....
        otherwise
            disp('Such a model has not yet been built in');
    end
    
    
    % Finally we wise to save 'tpm' into a format which is usable by the
    % IIT software.
    if print
        filename = sprintf('External Software/iit/tpm_matrix_M%d_K%.2f_%.2f_%.2f_P%.2f_%.2f.mat', model, k(1), k(2), k(3), p(1), p(2));
        save(filename, 'tpm');
    end
end

