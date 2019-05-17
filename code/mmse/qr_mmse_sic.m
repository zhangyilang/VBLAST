% function qr_mmse_sic
% description
%

function dec = qr_mmse_sic(r,H,ModType,sigma)
   
    [~,Nt] = size(H);
        
    dec = zeros(1,length(r));
    d   = zeros(1,Nt);
    z   = zeros(1,Nt);
    c   = zeros(1,Nt);
    
    Hext = [H;sigma.*eye(Nt)];      % extended format of Channel Matrix H 
    [Q,R] = qr(Hext);
    % [Q1,R1] = qr(Hext,0);
    r_ext = [r;zeros(Nt,1)];
    y = Q'*r_ext;
    
    for k = Nt:-1:1
        for i = (k+1):Nt
            d(k) = d(k) + R(k,i)*c(i);
        end
        z(k) = y(k) - d(k); 
        dec(k) = qamdemod(z(k)/R(k,k), ModType);
        c(k) = qammod(dec(k), ModType);
    end
end
