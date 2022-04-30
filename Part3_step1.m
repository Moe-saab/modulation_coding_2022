%% Author Mohamad Saab
%% Impact of CFO on bit error rate (only ISI)
% Initialization
nbits = 1000;
Nbps = 2;
cutoffreq = 2e9;
rolloff = 0.3;
m = 4; %upsampling ratio
Ntaps =100; %no of taps 
EbN0 = -4:16; %snr
Tsymb = 1/(2*cutoffreq);
symbrate = 1/Tsymb;
Fs = symbrate*m;	%sampling frequency
repeatnumber = 10;
ppm = [0 2 10 20 50 100];
CFO = ppm*cutoffreq*1e-6;
for repeat = 1:repeatnumber
    
    %bit generation
    bits_tx = randi(2,1,nbits)-1; %generate sequence of 0's and 1's
    
    %mapping
    if Nbps>1
	signal_tx = mapping(bits_tx.', Nbps, 'qam').';
    else
	signal_tx = mapping(bits_tx.', Nbps, 'pam').';
    end
    
    %upsampling
    upsampled_signal = upsample(signal_tx,m);
    
    %hRRC at Tx
    [hrrc_time,hrrc_frequency] = hrrc(Fs,Tsymb,Ntaps,rolloff);
    filtered_signal_tx = conv(upsampled_signal,hrrc_time);
    
    %adding noise
    signal_energy = trapz(abs(filtered_signal_tx).^2)*(1/Fs);
    Eb = signal_energy/(2*nbits);
    N0 = Eb./(10.^(EbN0/10));
    noise_power = 2*N0*Fs;
    
    adwg_noise = zeros(length(EbN0),length(signal_tx)*m+Ntaps-1);
    signal_rx = zeros(length(EbN0),length(signal_tx)*m+Ntaps-1);

    for j = 1:length(EbN0)
        adwg_noise(j,:) = sqrt(noise_power(j)/2).*(randn(1,length(signal_tx)*m+Ntaps-1)+1i*randn(1,length(signal_tx)*m+Ntaps-1));
        signal_rx(j,:) = filtered_signal_tx + adwg_noise(j,:);
    end
    
    
    % adding CFO (frequency offset)
    t1 = ((0:length(signal_rx)-1))*1/Fs;
    signal_rx_sync_errors=zeros(length(EbN0),size(signal_rx,2),length(CFO));
    
    for k=1:length(CFO)
        for i = 1:length(EbN0)
            signal_rx_sync_errors(i,:,k) = signal_rx(i,:).*exp(1j*(2*pi*CFO(k).*t1));
        end
    end
    
    %so you should not multiply with an exponential but with sine for
    %imaginary and cosine for the real part since we are working with
    %baseband
    % never mind this note
    % you can simulate the cfo error will cause a circle ring
    % increase a lot the sampling when simulating the timeshift
    % and on the receiver you don't sample 
    
    %hRRC at Rx 
    filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*m,length(CFO));
    for k = 1:length(CFO)
        for i =1:length(EbN0)
            filtered_signal_rx = conv(signal_rx_sync_errors(i,:,k),fliplr(hrrc_time));
            filtered_signal_rx(i,:,k) = filtered_signal_rx(Ntaps:end-(Ntaps-1));
        end                                                                                     %
    end 
    
end
