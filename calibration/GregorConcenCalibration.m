
% Gregor concentration calibration
APpos = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
GregorBcdConcen_nM = [53 40 23 14 12 6 4 2.5 1 0.5 0];  %Bcd gradient in nM

%1nM conversion to molecules per cubic micron (molecules per femtoliter)
oneNanoMolar = (10^(-6)) * (10^(-6))^3 * 6.022*(10^23); 

% Bcd gradient in molecules per cubic micron
GregorBcdConcen_molecPerCubeMicron = GregorBcdConcen_nM * oneNanoMolar

% 1 fL = 1 um^3 ~ excitation volume ~ 1 pixel
% We are integrating over 9 pixels to get the protein signal at the locus
numPixels = 9;
Gregor_totalBcdPerDisk = GregorBcdConcen_molecPerCubeMicron * numPixels


