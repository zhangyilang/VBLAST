% VBLAST Detection
%% Clear
clear
clc
close all

%% Parameters
Nt = 4;             % send antenna number
Nr = 4;             % receive antenna number

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
dec_mmse        = zeros(L,1);
dec_mmse_sorted = zeros(L,1);
dec_mmse_sqrd   = zeros(L,1);
   
EB_mmse         = zeros(1,length(SNRindB));
EB_mmse_sorted  = zeros(1,length(SNRindB));
EB_qr_mmse      = zeros(1,length(SNRindB));
EB_mmse_sqrd    = zeros(1,length(SNRindB));

for index = 1:length(SNRindB)
    SNR(index) = 10^(SNRindB(index)/10);
    sigma = sqrt(1/SNR(index));
    
    for simcnt= 1:SimTimes
        txMsgBits  = randi([0,1], [log2(ModType),L]);
        txMsgBitsInt = bi2de(txMsgBits', 'left-msb')';
        txMapped = qammod(txMsgBitsInt, ModType);
        
        x = reshape(txMapped,Nt,L/Nt);                                    % reshape data to Nt antennas
        AWGN_noise = sqrt(1/2)*sigma*(randn(Nr,L/Nt)+1j*randn(Nr,L/Nt));  % the Nois Matrix
        
        % ======== do detection procedure ========
        H = sqrt(1/Nt)*sqrt(1/2)*(randn(Nr,Nt) + 1j*randn(Nr,Nt));    % fast fading Rayleigh channel
        r = H*x + AWGN_noise;                                         % received signal
        for col_idx = 1:L/Nt

            rsic = r(:,col_idx);
            dec_mmse((col_idx-1)*Nt+1:col_idx*Nt) = vblast_mmse(rsic,H,ModType,sigma);
            dec_qr_mmse((col_idx-1)*Nt+1:col_idx*Nt) = qr_mmse_sic(rsic,H,ModType,sigma);
            dec_mmse_sqrd((col_idx-1)*Nt+1:col_idx*Nt) = qr_mmse_sic_sorted_v3(rsic,H,ModType,sigma);
            dec_mmse_sorted((col_idx-1)*Nt+1:col_idx*Nt) = vblast_mmse_sorted(rsic,H,ModType,sigma);
            
        end % end of dection loop
        det_mmse_bin = reshape(de2bi(dec_mmse,2,'left-msb')',1,[]);
        dec_qr_mmse_bin = reshape(de2bi(dec_qr_mmse,2,'left-msb')',1,[]);
        dec_mmse_sqrd_bin = reshape(de2bi(dec_mmse_sqrd,2,'left-msb')',1,[]);
        dec_mmse_sorted_bin = reshape(de2bi(dec_mmse_sorted,2,'left-msb')',1,[]);
        
        EB_mmse(index)    = EB_mmse(index) + sum(abs(det_mmse_bin~=txMsgBits(:)'));
        EB_qr_mmse(index) = EB_qr_mmse(index) + sum(abs(dec_qr_mmse_bin~=txMsgBits(:)'));
        EB_mmse_sqrd(index) = EB_mmse_sqrd(index) + sum(abs(dec_mmse_sqrd_bin~=txMsgBits(:)'));
        EB_mmse_sorted(index) = EB_mmse_sorted(index) +sum(abs(dec_mmse_sorted_bin~=txMsgBits(:)'));
    
    end %end of simcnt loop
end % end of SNR loop

TotalBits = ((L*log2(ModType))*SimTimes);
BER_mmse = EB_mmse./TotalBits;
BER_mmses_sorted = EB_mmse_sorted./TotalBits;
BER_qr_mmse = EB_qr_mmse./TotalBits;
BER_mmse_sqrd = EB_mmse_sqrd./TotalBits;

%% show the results. MMSE
figure(100);
semilogy(EbN0indB,BER_mmse,'-ro','LineWidth',2);hold on;
semilogy(EbN0indB,BER_qr_mmse,'-b*','LineWidth',2);hold on;
semilogy(EbN0indB,BER_mmse_sqrd,'-kv','LineWidth',2);hold on;
semilogy(EbN0indB,BER_mmses_sorted,'-mpentagram','LineWidth',2);hold on;
xlabel('Eb/N_{0} in dB');ylabel('BER');
legend('MMSE','MMSE-QRD','MMSE-SQRD','MMSE-BLAST');
grid on;
