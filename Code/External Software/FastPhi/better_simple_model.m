Prob = [0.8*0.75 0.2*0.75 0.25*0.8 0.2*0.25];
pmat = shuffle([.4 .3 .2]);


%%
%Mechanism which assigns always assigns the minimum amount of labor to a
%project. This mechanism will ignore a project instead of wasting a
%multiskilled team on a focused project

Mech1 = zeros([8 4]);
Mech1(:,1) = 0:7;

for i=0:7
    for P = 1:3
        Mech1(i+1, P+1) = Mech1(i+1,1);
        
        if ~bitget(i,3)
            if P==1 || (P==3 && bitget(i,1) && ~bitget(i,2))
                Mech1(i+1,P+1) = Mech1(i+1,P+1)+4;
            end
        end

        if ~bitget(i,2)
            if P==2 || (P==3 && bitget(i,1) &&~bitget(i,3))
                Mech1(i+1,P+1) = Mech1(i+1,P+1) + 2;
            end
        end

        if ~bitget(i,1)
            if P==3
                Mech1(i+1,P+1) = Mech1(i+1,P+1) + 1;
            end
        end
    end
end

tmat1 = zeros(8);

for i=1:8
    for j = 1:8
        for Pind=1:4
            ind = int8(Mech1(j,Pind))+1;
            tmat1(i,ind)=tmat1(i,ind)+pmat(i,j)*Prob(Pind);
        end
    end
end

disp(tmat1)



%%
%This mechanism will complete a project iff there is a perfect fit
%available. Thus it will not make a multiskilled team work on a focused
%project, and similarly it will not make multiple specialists work together
%on a multiskilled project.

Mech2 = zeros([8 4]);
Mech2(:,1) = 0:7;

for i=0:7
    for P = 1:3
        Mech2(i+1, P+1) = Mech2(i+1,1);
        
        if ~bitget(i,3)
            if P==1 
                Mech2(i+1,P+1) = Mech2(i+1,P+1)+4;
            end
        end

        if ~bitget(i,2)
            if P==2 
                Mech2(i+1,P+1) = Mech2(i+1,P+1) + 2;
            end
        end

        if ~bitget(i,1)
            if P==3
                Mech2(i+1,P+1) = Mech2(i+1,P+1) + 1;
            end
        end
    end
end




%disp(pmat);

tmat2 = zeros(8);

for i=1:8
    for j = 1:8
        for Pind=1:4
            ind = int8(Mech2(j,Pind))+1;
            tmat2(i,ind)=tmat2(i,ind)+pmat(i,j)*Prob(Pind);
        end
    end
end

disp(tmat2)



%%

%Mechanism which assigns always assigns the minimum amount of labor to a
%project. This mechanism will ignore a project instead of wasting a
%multiskilled team on a focused project

Mech3 = zeros([8 4]);
Mech3(:,1) = 0:7;

for i=0:7
    for P = 1:3
        Mech3(i+1, P+1) = Mech3(i+1,1);
        
        if ~bitget(i,3)
            if P==1 || (P==3 && bitget(i,1) && ~bitget(i,2))
                Mech3(i+1,P+1) = Mech3(i+1,P+1)+4;
            end
        end

        if ~bitget(i,2)
            if P==2 || (P==3 && bitget(i,1) &&~bitget(i,3))
                Mech3(i+1,P+1) = Mech3(i+1,P+1) + 2;
            end
        end

        if ~bitget(i,1)
            if P==3 || (bitget(i,2) && bitget(i,3))
                Mech3(i+1,P+1) = Mech3(i+1,P+1) + 1;
            end
        end
    end
end

tmat3 = zeros(8);

for i=1:8
    for j = 1:8
        for Pind=1:4
            ind = int8(Mech3(j,Pind))+1;
            tmat3(i,ind)=tmat3(i,ind)+pmat(i,j)*Prob(Pind);
        end
    end
end

disp(tmat3)




