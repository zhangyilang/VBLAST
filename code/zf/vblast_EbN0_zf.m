% VBLAST Detection
%% Clear
clear
clc
close all

%% Parameters
Nt = 4;             % send antenna number (set 1 to simulate SISO)
Nr = 4;             % receive antenna number (set 1 to simulate SISO)

L = 100;            % frame length
SimTimes = 1e3;     % static count per SNR (repeat times)

EbN0indB = 0:2:30;       % define Eb/N0 range
ModType = 4;             % modulation mode: 1=BPSK, 4=QPSK, 16=16QAM, 64=64QAM
% C/N = SNR = (Eb/No)*(fb/B) => SNR = 10log10(Eb/N0) + 10log10(fb/B)
SNRindB = EbN0indB + 10*log10(log2(ModType));

% Average Symbol Energy 1 BPSK, 2 QPSK, 5 16QAM, etc...
SNR = zeros(1,length(SNRindB));


%% Loop
% initialization
dec_zf          = zeros(L,1);
dec_zf_sorted   = zeros(L,1);
dec_qr_zf       = zeros(L,1);
dec_zf_sqrd     = zeros(L,1);
   
EB_zf           = zeros(1,length(SNRindB));
EB_zf_sorted    = zeros(1,length(SNRindB));
EB_qr_zf        = zeros(1,length(SNRindB));
EB_zf_sqrd      = zeros(1,length(SNRindB));

for index = 1:length(SNRindB)
    SNR(index) = 10^(SNRindB(index)/10);
    sigma = sqrt(1/SNR(index));
    
    for simcnt= 1:SimTimes
        txMsgBits  = randi([0,1], [log2(ModType),L]);
        txMsgBitsInt = bi2de(txMsgBits', 'left-msb')';
        txMapped = qammod(txMsgBitsInt, ModType);
        
        x = reshape(txMapped,Nt,L/Nt);                                    % reshape data to Nt antennas
        % x = sqrt(1/Es)*x;                                               % normalization 
        AWGN_noise = sqrt(1/2)*sigma*(randn(Nr,L/Nt)+1j*randn(Nr,L/Nt));  % Gaussian noise Matrix
        
        % ======== do detection procedure ========
        H = sqrt(1/Nt)*sqrt(1/2)*(randn(Nr,Nt) + 1j*randn(Nr,Nt));    % fast fading Rayleigh channel
        r = H*x + AWGN_noise;                                         % received signal
        for col_idx = 1:L/Nt

            rsic = r(:,col_idx);
            dec_zf((col_idx-1)*Nt+1:col_idx*Nt)   = vblast_zf(rsic,H,ModType);
            dec_qr_zf((col_idx-1)*Nt+1:col_idx*Nt) = qr_zf_sic(rsic,H,ModType);
            dec_zf_sorted((col_idx-1)*Nt+1:col_idx*Nt) = vblast_zf_sorted(rsic,H,ModType);
            dec_zf_sqrd((col_idx-1)*Nt+1:col_idx*Nt) = qr_zf_sic_sorted(rsic,H,ModType);
            
        end
         dec_zf_bin = reshape(de2bi(dec_zf,2,'left-msb')',1,[]);
         dec_zf_sorted_bin = reshape(de2bi(dec_zf_sorted,2,'left-msb')',1,[]);
         dec_qr_zf_bin = reshape(de2bi(dec_qr_zf,2,'left-msb')',1,[]);
         dec_zf_sqrd_bin = reshape(de2bi(dec_zf_sqrd,2,'left-msb')',1,[]);
        
         EB_zf(index)      = EB_zf(index) + sum(abs(dec_zf_bin~=txMsgBits(:)'));
         EB_zf_sorted(index)   = EB_zf_sorted(index) + sum(abs(dec_zf_sorted_bin~=txMsgBits(:)'));      
         EB_qr_zf(index)   = EB_qr_zf(index) + sum(abs(dec_qr_zf_bin~=txMsgBits(:)'));
         EB_zf_sqrd(index) = EB_zf_sqrd(index) + sum(abs(dec_zf_sqrd_bin~=txMsgBits(:)'));
        
    end %end of simcnt loop
end % end of SNR loop

TotalBits = ((L*log2(ModType))*SimTimes);
 BER_zf   = EB_zf./TotalBits;
 BER_zf_sorted  = EB_zf_sorted./TotalBits;
 BER_qr_zf = EB_qr_zf./TotalBits;
 BER_zf_sqrd = EB_zf_sqrd./TotalBits;

%% show results
figure(1);
semilogy(EbN0indB,BER_zf,'-ro','LineWidth',2);hold on;
semilogy(EbN0indB,BER_qr_zf,'-b*','LineWidth',2);hold on;
semilogy(EbN0indB,BER_zf_sqrd,'-kv','LineWidth',2);hold on;
semilogy(EbN0indB,BER_zf_sorted,'-mpentagram','LineWidth',2);hold on;
xlabel('E_b/N_{0} in dB');ylabel('BER');
legend('ZF','ZF-QRD','ZF-SQRD','ZF-BLAST');
grid on;
