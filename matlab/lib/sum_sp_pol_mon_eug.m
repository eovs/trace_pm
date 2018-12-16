function  [wx,dx,cx,xs]=sum_sp_pol_mon(wc,dc,cc,cs,   wb,db,cb,bs)
% x,y are sparse polynomials in form (weight,degrees, coefs.
% z =x+y; 
if wc==0 
    wx=wb; dx=db; cx=cb; xs=bs;
    return; 
end
if wb==0
    wx=wc; dx=dc; cx=cc; xs=cs;
    return; 
end

mywx = max( dc(wc), db(wb));

polc = zeros( 1, mywx );
polb = zeros( 1, mywx );
symc = zeros( 1, mywx );
symb = zeros( 1, mywx );
 
if dc(1) == 0 
    i = 0;
end

% mywx
% dc
% cc
polc(dc(1:wc)+1) = cc(1:wc);
polb(db(1:wb)+1) = cb(1:wb);
symc(dc(1:wc)+1) = cs(1:wc);
symb(db(1:wb)+1) = bs(1:wb);
% 
% polr = polc + polb;
% symr = symc + symb;

% polt = polc .* polb;
% f = find( polt ~= 0 )
% mm = min( symc, symb );
% 
% symr(f) = mm(f);
% mydx = find( polr ~= 0  );
% mycx = polr( dx );
% myxs = symr( dx );
%
%mywx = sum( polr ~= 0 );


% artificial extra degrees for termination 
dc(wc+1)=db(wb)+1;  
db(wb+1)=dc(wc)+1;  


% Merge 
wx=wb+wc;
dx=zeros(1,wx);
cx=zeros(1,wx);
xs=zeros(1,wx);

ic=0; ib=0; % start with smallest degrees
wx=0; 
while ic<wc || ib<wb
    if db(ib+1)<dc(ic+1)
        % read B and increase index
        wx=wx+1; ib=ib+1; 
        dx(wx)=db(ib);
        cx(wx)=cb(ib);
        xs(wx)=bs(ib);
    else
        if db(ib+1)>dc(ic+1)
            wx=wx+1;ic=ic+1; 
            dx(wx)=dc(ic);
            cx(wx)=cc(ic);
            xs(wx)=cs(ic);
        else  % db==dc
            a=cb(ib+1)+cc(ic+1);
            ic=ic+1; ib=ib+1; 
            if a~=0
                wx=wx+1;           
                cx(wx)=a;
                dx(wx)=db(ib);
                xs(wx)=min(cs(ic), bs(ib));
            end
        end
    end
end
% remove extra zeros if necessary
dx=dx(1:wx);
cx=cx(1:wx);
xs=xs(1:wx);

if wx ~= mywx
    mywx = 0;
end;