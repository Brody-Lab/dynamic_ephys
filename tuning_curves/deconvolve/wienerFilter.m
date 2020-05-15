function [yhat H] = wienerFilter(ideal,observation);
R=1;     % response function     % ... no smearing

% estimate noise from ideal
noise = observation-ideal;
% work out how long to make FFT
N=length(observation);

% Wiener filter
Sf2=real(fft(ideal,N*2-1)).^2;   % Smeared ideal
Nf2=real(fft(noise,N*2-1)).^2;   % noise
Cf=real(fft(observation,N*2-1)); % ~= sqrt(Sf2+Nf2); % Corrupted ideal
H=Sf2./(Sf2+Nf2);              % Optimal filter
Yhat=H.*Cf/R;                  % 'uncorrupted' ideal estimate ...
yhat=real(ifft(Yhat));           % ..... in time domain
% ...compensate for FFT being two sided in matlab   
yhat=yhat(1:length(observation)); 


