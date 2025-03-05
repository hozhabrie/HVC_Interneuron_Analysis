file = abfload('C:\Users\User\Dropbox (NYU Langone Health)\Int Juxta\eh12\c4_0\c4_0_filteredspikes.ABF');

fs = 50000;
spikes = file(:,1); 
axis(1) = subplot(2,1,1)
plot(linspace(0, length(file)/fs, length(file)),spikes)
axis(2) = subplot(2,1,2)

%%
clc; close all;

% figure(1)
% s1 = spikes(50*fs:60*fs-1);
% plot(s1)
% shg

s1 = spikes(50*fs:60*fs-1);
s2 = spikes(70*fs:80*fs-1);
s3 = spikes(135*fs:145*fs-1);
s4 = spikes(182*fs:192*fs-1);
s5 = spikes(208*fs:238*fs-1);

allSpont = s5;%[s1, s2, s3, s4, s5];

for i = 1%:size(allSpont,2)
    plot(linspace(0,length(allSpont)/fs, length(allSpont)), allSpont(:,i)-2.5*i, 'k', 'LineWidth', .5);
    hold on
end

