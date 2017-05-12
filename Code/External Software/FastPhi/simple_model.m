T = zeros(8);

p = 0.2;
q = 0.25;

k1 = 0.2;
k2 = 0.3;
k3 = 0.4;

T(1,1)=(1-p)*(1-q);
T(1,2) = p*q;
T(1,3) = q*(1-p);
T(1,4) = 0;
T(1,5) = p*(1-q);
T(1,6) = 0;
T(1,7) = 0;
T(1,8) = 0;

check = sum(T(1,:));
disp(check);

T(2,1) = k3*(1-p)*(1-q);
T(2,2) = (1-k3)*(1-p)*(1-q) + k3*p*q;
T(2,3) = k3*q*(1-p);
T(2,4) = (1-k3)*q*(1-p);
T(2,5) = k3*p*(1-q);
T(2,6) = (1-k3)*p*(1-q);
T(2,7) = 0;
T(2,8) = (1-k3)*p*q;

check = sum(T(2,:));
disp(check);

row = 3;

T(row,1) = k2*(1-p)*(1-q);
T(row,2) = k2*p*q;
%The third term gives rise to waste
T(row,3) = (1-k2)*(1-p)*(1-q) + k2*q*(1-p) + (1-k2)*q*(1-p);
T(row,4) = (1-k2)*p*q;
T(row,5) = k2*p*(1-q);
T(row,6) = 0;
T(row,7) = (1-k2)*p*(1-q);
T(row,8) = 0;

check = sum(T(row,:));
disp(check);

%%
row = 4;

T(row,1) = k2*k3*(1-p)*(1-q);
%waste
T(row,2) = k2*(1-k3)*(1-p)*(1-q) + k2*k3*p*q + k2*(1-k3)*p*q;
%waste
T(row,3) = (1-k2)*k3*(1-p)*(1-q) + k2*k3*q*(1-p) + (1-k2)*k3*q*(1-p);
T(row,4) = (1-k2)*(1-k3)*(1-p) + k2*(1-k3)*q*(1-p) + (1-k2)*k3*p*q;
T(row,5) = k2*k3*p*(1-q);
T(row,6) = k2*(1-k3)*p*(1-q);
T(row,7) = k3*(1-k2)*p*(1-q) ;
T(row,8) = (1-k2)*(1-k3)*p;

check = sum(T(row,:));
disp(check);

%%
row = 5;

T(row,1) = k1*(1-p)*(1-q);
T(row,2) = k1*p*q;
T(row,3) = k1*q*(1-p);
T(row,4) = 0;
T(row,5) = k1*p*(1-q) + (1-k1)*(1-q)*(1-p);
T(row,6) = (1-k1)*p*q;
T(row,7) = (1-k1)*q*(1-p) + (1-k1)*p*(1-q);
T(row,8) = 0;

check = sum(T(row,:));
disp(check);

%%

row = 6;

T(row,1) =
T(row,2) = 
T(row,3) = 
T(row,4) = 
T(row,5) = 
T(row,6) = 
T(row,7) = 
T(row,8) = 

check = sum(T(row,:));
disp(check);


row = 7;

T(row,1) =
T(row,2) = 
T(row,3) = 
T(row,4) = 
T(row,5) = 
T(row,6) = 
T(row,7) = 
T(row,8) = 

check = sum(T(row,:));
disp(check);


row = 8;

T(row,1) =
T(row,2) = 
T(row,3) = 
T(row,4) = 
T(row,5) = 
T(row,6) = 
T(row,7) = 
T(row,8) = 

check = sum(T(row,:));
disp(check);


