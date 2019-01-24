%%
clear
clc

% Article: "Principal eigenvector localization and centrality in networks: revisited"
% Priodyuti Pradhan, Angeliya C.U. and Sarika jalan.


%Algorithm details is in "Eigenvalue crossing in principal eigenvector
% localized networks", Priodyuti Pradhan and Sarika Jalan, arXiv:1810.00243v1, 2018.

% email: priodyutipradhan@gmail.com, sarikajalan9@gmail.com
 
%N = 10000; E = 16000; 
%N = 10922; E = 18061; 
%N = 10498; E = 52490; 
%N = 20422; E = 163376; 
%N = 1133;
%k = 9.62
%E = ceil(k*N)/2

%%
%Input: Given a value N and E such that 
% E <(N-1)+ceil(((N+2)*sqrt(3*(N+2)))/9)

%Output: For k1, n1 and n2. Similarly
%        For k2, n1 and n2.

N = 2408
E = 14806

%N = 32768
%E = 1391268

%E = (N + ((N+2)*sqrt(3*(N+2)))/9)

if E> (N-1)+ceil(((N+2)*sqrt(3*(N+2)))/9)
 fprintf('Discriminant is greater than zero!!! Can not use roots!\n')
 return
end

eps2 = 0.00002;
sigma = 1.0 - eps2;

p = -2*sigma - 4;
q = 8*sigma + sigma^2 + 1 - N;
r = 2*E - 4*sigma^2;


a = (1/3) * (3*q - p*p);
b = (1/27) * (2*p*p*p - 9*p*q + 27*r);


T1 = b*b/4+a*a*a/27;

if T1==0.0 || T1>0.0
   A = nthroot(((-b/2) + sqrt(T1)),3);
   B = nthroot(((-b/2) - sqrt(T1)),3);
else
   A = ((-b/2) + sqrt(T1))^(1/3);
   B = ((-b/2) - sqrt(T1))^(1/3);
end

y1 = A + B;
y2 = -(1/2)*(A+B) - (i*sqrt(3)/2)*(A-B);
y3 = -(1/2)*(A+B) + (i*sqrt(3)/2)*(A-B);

%%y2 = -(1/2)*(A+B) + ((sqrt(-1)*sqrt(3))/2)*(A-B);
%%y3 = -(1/2)*(A+B) - ((sqrt(-1)*sqrt(3))/2)*(A-B);


k1 = y1 - p/3;
k2 = y2 - p/3;
k3 = y3 - p/3;
 
fprintf('\n')
fprintf('x1: %0.2f\nx2: %0.2f\nx3: %0.2f\n',k1,k2,k3);

n1 = ceil((ceil(k1)+eps2-1)^2);
n2 = N-n1-1;
fprintf('n1: %d n2: %d\n',n1,n2);

n1 = ceil((ceil(k2)+eps2-1)^2);
n2 = N-n1-1;
fprintf('n1: %d n2: %d\n',n1,n2);

n1 = ceil((ceil(k3)+eps2-1)^2);
n2 = N-n1-1;
fprintf('n1: %d n2: %d\n',n1,n2);


