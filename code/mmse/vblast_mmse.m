% function vblast_mmse.m
% description :
%

function  dec = vblast_mmse(r,H,ModType,sigma)
    [~,Nt] = size(H);
    G = (H'*H+sigma.^2*eye(Nt))\H';
    y = G*r;
    dec = qamdemod(y, ModType);
    dec = dec';
end
% Matrix G : Nt*Nr

