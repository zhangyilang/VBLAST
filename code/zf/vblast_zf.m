% function vblast_zf.m
% description : lineral zero forcing 
%

function  dec = vblast_zf(r,H,ModType)
   
    G = (H'*H)\H';      % G = pinv(H)? Moore-Penrose pesudoinverse
    y = G*r;
    dec = qamdemod(y, ModType);
    dec = dec';
end
% Matrix G : Nt*Nr
