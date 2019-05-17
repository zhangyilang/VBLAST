% Implement MIMO VBLAST to transmit pictures
%% Clear
clear
clc
close all

%% Parameters
Nt = 4;             % send antenna number
Nr = 4;             % receive antenna number

EbN0indB = 10:5:20;     % Eb/N0 range
ModType = 4;            % modulation mode: 1=BPSK, 4=QPSK, 16=16QAM, 64=64QAM
% C/N = SNR = (Eb/No)*(fb/B) => SNR = 10log10(Eb/N0) + 10log10(fb/B)
SNRindB = EbN0indB + 10*log10(log2(ModType));

%% Read picture
img = imread('test.png');
pixels = length(img(:));
img_bin = de2bi(reshape(img, pixels, 1));        % convert into binary form
L = length(img_bin(:));         % frame length
txMsgBits = reshape(img_bin, log2(ModType), []);
txMapped = qammod(txMsgBits, ModType);
x = reshape(txMapped,Nt,L/Nt);  % reshape data to Nt antennas

%% Transmit
dec_zf_sqrd = zeros(L,length(SNRindB));
% waitbar
h = waitbar(0,'Evaluation In Process  (~£þ¨Œ£þ)~');
for i = 1:length(SNRindB)
    SNR = 10^(SNRindB(i)/10);
    sigma = sqrt(1/SNR);
    AWGN_noise = sqrt(1/2)*sigma*(randn(Nr,L/Nt)+1j*randn(Nr,L/Nt));  % noise matrix

    % ======== do detection procedure ========
    H = sqrt(1/Nt)*sqrt(1/2)*(randn(Nr,Nt) + 1j*randn(Nr,Nt));  % fast fading Rayleigh channel
    r = H * x + AWGN_noise;                                     % received signal
    for col_idx = 1:L/Nt
        dec_zf_sqrd((col_idx-1)*Nt+1:col_idx*Nt, i) = qr_zf_sic_sorted(r(:,col_idx),H,ModType);
        % update waitbar
        waitbar((i-1)/length(SNRindB)+col_idx*Nt/L/length(SNRindB));
    end
end
% close waitbar
close(h)

%% Receive and reconstruct
dec_zf_sqrd = mod(dec_zf_sqrd, 2);  % binary domain
img_r = cell(3, 1);                 % cell for storing recovered images
for i = 1:length(SNRindB)
    img_rec_bin = reshape(dec_zf_sqrd(:,i), size(img_bin));
    img_rec = bi2de(img_rec_bin, 2);
    img_rec = uint8(reshape(img_rec, size(img)));
    img_r{i} = img_rec;
end

%% show results
figure(1)
subplot(2,2,1)
imshow(img)
title('origin image')
subplot(2,2,2)
imshow(img_r{1})
title('SNR=10dB')
subplot(2,2,3)
imshow(img_r{2})
title('SNR=15dB')
subplot(2,2,4)
imshow(img_r{3})
title('SNR=20dB')
