%% LAND: Controller %%
%---------------------------------------------------------------------------------------------------------------------------------
clear
% close all
clc

Fs = 100;               % sampling frequency [Hz]
t = 0:(1/Fs):5;         % time vector [s]
n = length(t);          % # data points

f = linspace(0.5,3,n);  % frequency range [Hz]
I = sin(2*pi*f.*t);     % input signal (chirp)
delay = 0.1;            % delay [s]
d = delay*Fs;           % delay [# idx]
I(1:d) = 0;             % set input to zero for time of delay

PSM = [0.050 , 0.000 , 0.167 , 0.167 , 0.000 , 0.000 , 0.000 , 0.000 ];	% smooth propotional gains
VSM = [0.000 , 0.500 , 0.000 , 0.500 , 0.000 , 0.000 , 0.500 , 0.500 ];	% smooth velocity gains
PSA = [0.000 , 0.000 , 0.000 , 0.000 , 1.000 , 1.000 , 1.000 , 1.000 ];	% saccadic propotional gains
VSA = [0.000 , 0.000 , 0.000 , 0.000 , 0.000 , 5.000 , 2.000 , 2.000 ];	% saccadic velocity gains
LIM = [0     , 0     , 0     , 0     , 0     , 0     , 0     , 1     ];	% limit velocity (not used yet)

O = zeros(length(PSM),n); % preallocate output vector
E = zeros(length(PSM),n); % preallocate error vector

figure (1) ; clf ; hold on
for j = 1:length(PSM)
    for k = (d+2):n
        if (PSM(j) || VSM(j)) && (PSA(j) || VSA(j))
            O(j,k) = O(j,k-1) + PSM(j)*(I(k-d) - O(j,k-d)) + VSM(j)*((I(k-d) - O(j,k-d))-(I(k-1-d) - O(j,k-1-d)));       
            if mod(k,d)==0 && k<(n-d)
                O(j,k:(k+d)) = PSA(j)*I(k-d) + VSA(j)*(I(k-d) - I(k-d-1));
            end
        elseif PSM(j) || VSM(j)
            O(j,k) = O(j,k-1) + PSM(j)*(I(k-d) - O(j,k-d)) + VSM(j)*((I(k-d) - O(j,k-d))-(I(k-1-d) - O(j,k-1-d)));       
        elseif PSA(j) || VSA(j)
            if mod(k,d)==0 && k<(n-d)
                O(j,k:(k+d)) = PSA(j)*I(k-d) + VSA(j)*(I(k-d) - I(k-d-1));
            end
        end
    end
    O_out(j,:) = [O(j,d:n) , zeros(1,d-1)];
  	E(j,:) = I - O_out(j,:);

    h = figure (1);
    subplot(ceil(length(PSM)/2),2,j) ; hold on ;
        title({['PSM = ' num2str(PSM(j)) '  ' 'VSM = ' num2str(VSM(j))],...
               ['PSA = ' num2str(PSA(j)) '  ' 'VSA = ' num2str(VSA(j))]})
        plot(t,I)
        plot(t,O_out(j,:))
        plot(t,E(j,:))
        xlim([0 max(t)])
        box on
        if j==1 ; legend('Input','Output','Error') ; end       
end
