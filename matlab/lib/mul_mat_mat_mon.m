function  [WX,DX,CX,XS,LX]=mul_mat_mat_mon(WB,DB,CB,BS,A,AS,M)
% A,B are degree matrices
% B is 3D
% A is 2D
% M is modulo

[bb,~,LB]=size(DB);

LX=LB+1;
db=zeros(1,LB); % Column degrees 
cb=zeros(1,LB); % coefs
bs=zeros(1,LB); % ace

% ARRAY for result
WX=zeros(bb,bb);    % weights
DX=zeros(bb,bb,LX); % degrees
CX=zeros(bb,bb,LX); % coefs
XS=zeros(bb,bb,LX); % ACEs

% main loop
w=0;  % max weight
for I=1:bb % rows of B
    for J=1:bb % columns of A (square marices)
        wx=0; dx=0; cx=0; xs=0;      % accumulator
        for j=1:bb % loop along column
            if WB(I,j)>0 && A(j,J)>=0
            a=A(j,J); as=AS(j,J);
            wb=WB(I,j);
            for h=1:wb  % read from 3D arrays
                db(h)=DB(I,j,h);
                cb(h)=CB(I,j,h);
                bs(h)=BS(I,j,h);
            end
            % Multiplication 
            bs(1:wb)=bs(1:wb)+as;         % new col weights
            
            db(1:wb)=db(1:wb)+a;   % new degrees
                       
            t=find(db>=M,1);
             
             if ~isempty(t) 
             if t>1                 % Split
                 db1=db(1:t-1); db2=db(t:wb)-M;
                 cb1=cb(1:t-1); cb2=cb(t:wb);
                 bs1=bs(1:t-1); bs2=bs(t:wb);
                 db1(t)=M;      db2(wb-t+2)=M;  
                % Merge
                i1=0; i2=0; i=0;
                while i1<t-1 || i2<wb-t+1
                    if db1(i1+1)<db2(i2+1)
                        i=i+1; i1=i1+1;
                        db(i)=db1(i1);
                        cb(i)=cb1(i1);
                        bs(i)=bs1(i1);
                    else
                        if db1(i1+1)>db2(i2+1)
                            i=i+1; i2=i2+1;
                            db(i)=db2(i2);
                            cb(i)=cb2(i2);
                            bs(i)=bs2(i2);
                        else
                            aa=cb1(i1+1)+cb2(i2+1);
                            i1=i1+1; i2=i2+1;
                            if aa~=0
                                i=i+1;
                                cb(i)=aa;
                                db(i)=db1(i1);
                                bs(i)=min( bs1(i1), bs2(i2) );
                            end
                        end
                    end
                        
                        
                end
             else % t==1
                db=mod(db,M);
             end
             end
            % Accumulation
        
            [wx,dx,cx,xs]=sum_sp_pol_mon(wx,dx,cx,xs, wb,db,cb,bs);
            
            end
         end
         % write c to x
         WX(I,J)=wx;    % weights
         for h=1:wx
             DX(I,J,h)=dx(h); % degrees
             CX(I,J,h)=cx(h); % coefs
             XS(I,J,h)=xs(h); % ACEs
         end
         w=max(w,wx);  % max weight
    end
end
% Cut extra zeros
if w<max(LX,2)
   LX=w;
   DX=DX(:,:,1:w);
   CX=CX(:,:,1:w);
   XS=XS(:,:,1:w);
end
    
