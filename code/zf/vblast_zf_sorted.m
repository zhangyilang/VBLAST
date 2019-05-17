% function vblast_zf_sorted.m
% description : the ZF V-BLAST algorithm.
%

function  dec = vblast_zf_sorted(rsic,H,ModType)
%     thisMode   = modem.qammod(ModType);
%     thisDemod  = modem.qamdemod(ModType);
%     Es = (mean(thisMode.Constellation .* conj(thisDemod.Constellation))) ;
    
    dec = zeros(1,length(rsic));

    [~,Nt] = size(H);
    k = zeros(1,Nt);

    for i = 1:Nt
        if(i==1)
            G = pinv(H);
            [~,k0]=min(sum(abs(G).^2,2));
            r = rsic;
        end
        
        k(i) = k0;
        y = G(k(i),:)*r;
        dec(k(i)) = qamdemod(y, ModType);
        a = qammod(dec(k(i)), ModType);
        r = r - a*H(:,k(i));

        H(:,k(i)) = 0;
        G = pinv(H);
        temp = sum(abs(G).^2,2);
        temp(k(1:i)) = Inf;
        [~,k0] = min(temp);
    end
end
% Matrix G : Nt*Nr
% Vector w : Nr*1
% Vector r : Nr*1
