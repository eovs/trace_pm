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

