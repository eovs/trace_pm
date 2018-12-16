function [S,SA]=trace_bound_pol_mon_pm(HD,M,gmax)
% HP is polynomial biadjacency matrix of Tanner graph
% A is edge incidence poly-matrix (labeled)
% M is modulo
% N is number of edges (nonzero coefs in HD)


S=zeros(1,gmax); % Spectrum 
SA=zeros(1,gmax); % Spectrum 
[HP,HA]=tanner_mon(HD);
[A,AC]=hp2a_mon_pm(HP,HA,M);

N=size(HP,2);
WB=double(A>=0); % B  is accumulator WB are Hamming weights, 2D, integer 
DB=zeros(N,N,2);          % degrees,3D
DB(:,:,1)=A;                  % degrees, 3D
CB=zeros(N,N,2);          % coefficients, 3D
CB(:,:,1)=WB;                 %  coefs, 3D
BS=zeros(N,N,2); 
BS(:,:,1)=AC;                 % Accumullator for ACE, 3D
for d=4:2:gmax
    [WB,DB,CB,BS]=mul_mat_mat_mon(WB,DB,CB,BS,A,AC,M);
    for I=1:N
        wb=WB(I,I);
        if wb>0 && DB(I,I,1)==0
            % Increase spectra
            S(d)=S(d)+CB(I,I,1);
            if SA(d)>0
                SA(d)=min(SA(d),BS(I,I,1));
            else
                SA(d)=BS(I,I,1);
            end
            % Subtract constant
            wb=wb-1;
            WB(I,I)=wb;
            DB(I,I,1:wb)=DB(I,I,2:wb+1);
            CB(I,I,1:wb)=CB(I,I,2:wb+1);            
            %BS(I,I,1:wb)=BS(I,I,2:wb+1);
        end
    end
end

% 24, 24, 0, 48, 72, 264

for i=4:2:gmax
    S(i)=S(i)/i;
end
SA=SA/2;

