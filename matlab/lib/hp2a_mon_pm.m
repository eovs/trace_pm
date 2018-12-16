%function [A,WA,Am,WAm,Ap,WAp,Apm,WApm]=hp2a_mon_pm(HD,W,M)
function [Apm,WApm]=hp2a_mon_pm(HD,W,M)
% Construct polynomial edge-adjacency matrix A
% from Tanner graph
% HP is a polynomial parity-check matrix

% M is module

N=size(HD,2); 
B=HD>=0;  % base graph
A=-ones(2*N,2*N); 
WA=zeros(2*N,2*N); 
% minus
Am=-ones(N,N); 
WAm=zeros(N,N); 
% plus
Ap=-ones(N,N); 
WAp=zeros(N,N); 


% Edges: loop over columns of Tanner
E=zeros(2*N,2);
X=zeros(2*N,1);  % degrees
S=zeros(2*N,1);  % column weights
for i=1:N
    f=find(B(:,i)>0);
    E(i,1:2)=f(1:2);   % forward  part 1 to part 2
    E(i+N,1:2)=f([2 1]); % back  from 2 to 1;
    X(i)=HD(f(2),i);
    X(i+N)=mod(M-X(i),M);
    S(i)=W(f(2),i);
    S(i+N)=W(f(2),i);
end
% Edge adjacency matrix (Aminus)
for i=1:N          % i is head of the current edge
    j=E(i,2);        % j is tail of the current edge
    % find all edges starting at j ending anywhere except i+N or i-N 
    f=find(E(N+1:2*N,1)==j);
    f=f+N;
    for h=1:length(f)
        if f(h)~=(i+N)
             A(i,f(h))=X(f(h));
             Am(i,f(h)-N)=X(f(h));
             WA(i,f(h))=S(f(h));
             WAm(i,f(h)-N)=S(f(h));
        end
    end
end

%    (Aplus)
for i=N+1:2*N          % i is head of the current edge
    j=E(i,2);        % j is tail of the current edge
    % find all edges starting at j ending anywhere except i+N or i-N 
    % if i<=N, ib=i+N; else ib=i-N; end
    f=find(E(1:N,1)==j); 
    for h=1:length(f)
        if f(h)~=i-N
             A(i,f(h))=X(f(h));
             Ap(i-N,f(h))=X(f(h));
             WA(i,f(h))=S(f(h));
             WAp(i-N,f(h))=S(f(h));
        end
    end
end


% Multiplication Am*Ap
Apm=-ones(N);
WApm=zeros(N);
for I=1:N   % rows Am
    for J=1:N % columns Ap
        %s=-1; ws=0;  % accumulator
        for j=1:N % along row
            if Am(I,j)>=0 && Ap(j,J)>=0
                Apm(I,J)=mod(Am(I,j)+Ap(j,J),M); 
                WApm(I,J)=WAm(I,j)+WAp(j,J);
            end
        end
    end
end

    
     
    





return;

