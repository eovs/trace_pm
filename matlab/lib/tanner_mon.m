function [T,W]=tanner_mon(H)
% Generate  branch labels for the Tanner 
% graph of the code with parity check matrix H

[r,n]=size(H);
B=H>=0;   % base matrix
cw=sum(B);  % column weights
T=-ones(n+r, sum(sum(B)));
W=zeros(n+r, sum(sum(B)));
IC=0; % current column
for i=1:n
    f=find(B(:,i));
    for j=1:cw(i)
        IC=IC+1;
        T(i,IC)=1;
        T(n+f(j),IC)=H(f(j),i);
        W(n+f(j),IC)=cw(i);
    end;
end;
