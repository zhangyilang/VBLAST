% function vblast_mmse_sorted.m
% description :
%

function  dec = vblast_mmse_sorted(rsic,H,ModType,sigma)
    
    dec = zeros(1,length(rsic));

    [~,Nt] = size(H);
    k = zeros(1,Nt);

    for i = 1:Nt
        % do initialzation
        if(i==1)
            G = (H'*H+sigma.^2*eye(Nt))\H';
            [~,k0]=min(sum(abs(G).^2,2));
            r = rsic;
        end
        
        k(i) = k0;
        y = G(k(i),:)*r;                % 这个地方作�?�的写法有不�?,作�?�公式矩阵相乘，维度对不�?,这里多做了一次向量转�?
        dec(k(i)) = qamdemod(y, ModType);
        a = qammod(dec(k(i)), ModType);
        r = r- a*H(:,k(i));

        H(:,k(i)) = 0;
        G = (H'*H+sigma.^2*eye(Nt))\H';
        temp = sum(abs(G).^2,2);
        temp(k(1:i)) = Inf;
        [~,k0] = min(temp);
    end
end
% Matrix G : Nt*Nr
% Vector w : Nr*1
% Vector r : Nr*1

