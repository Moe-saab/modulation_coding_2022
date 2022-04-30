%part1
nbits = 1000;
Nbps = 2;
cutoffreq = 1e6;
rolloff = 0.3;
m = 4; %upsampling ratio
Ntaps =100; %no of taps 
EbN0 = 6; %snr
Tsymb = 1/(2*cutoffreq);
symbrate = 1/Tsymb;
Fs = symbrate*m;	%sampling frequency


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
% [h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,rolloff,Nbps,1,m);
[hrrc_time,hrrc_frequency] = hrrc(Fs,Tsymb,Ntaps,rolloff);
filtered_signal_tx = conv(upsampled_signal,hrrc_time);
figure("Name","TX signal");
subplot(1,3,1);
plot(signal_tx,'g*');
title("Symblos to be transmitted before upsampling");
subplot(1,3,2);
plot(upsampled_signal,'b*'); 
title("upsampled modulated signal","before filtering")
subplot(1,3,3);
plot(filtered_signal_tx,'r*'); title("filtered signal to be transmited");

%adding noise
signal_energy = trapz((abs(filtered_signal_tx).^2))*(1/Fs);
Eb = signal_energy/(2*nbits);
N0 = Eb./(10.^(EbN0/10));
noise_power = 2*N0*Fs;
% num of symbols multiplied by the upsampling rate + number of taps
adwg_noise = sqrt(noise_power/2).*(randn(1,nbits/Nbps * m +Ntaps-1) +1i*randn(1,nbits/Nbps * m +Ntaps-1));
signal_rx = filtered_signal_tx + adwg_noise;

%hRRC at Rx 
filtered_signal_rx = conv(signal_rx, fliplr(hrrc_time));
filtered_signal_rx = filtered_signal_rx(1,Ntaps:end-(Ntaps-1)); %didn't know why

%downsampling
downsampled_signal = downsample(filtered_signal_rx, m);

%figures
figure("Name","RX signal");
subplot(1,3,1);
plot(signal_rx,"r*")
title("received signal before filtering");
subplot(1,3,2);
plot(filtered_signal_rx,"y*")
title("Filtered RX signal");
subplot(1,3,3);
plot(downsampled_signal,"g*");
title("Downsampled RX signal");


%demapping
if Nbps>1
	bits_rx = demapping(downsampled_signal.', Nbps, 'qam');
else
	bits_rx = demapping(real(downsampled_signal.'), Nbps, 'pam');
end


%BER calculation
BER = 0;
for i=1:nbits
	if(bits_tx(i)~=bits_rx(i))
		BER = BER + 1;
	end
end
disp("the error bits are")
disp(BER);
disp("So BER in percent")
disp(BER/nbits)