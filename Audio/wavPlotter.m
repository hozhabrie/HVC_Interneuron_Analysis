clear all; close all; clc;
baseFolder='V:\Ellie\gcamp\dlx26_trainingrig\20190110\';
fileNames = dir([baseFolder '*.wav']);

for i = 1:length(fileNames)
    fname = [baseFolder fileNames(i).name];
    [audioFileRaw, fs] = audioread(fname);
    audioFileNorm = audioFileRaw./max(abs(audioFileRaw));
    filename=fileNames(i).name;
    figure;
    vigiSpec(audioFileNorm,fs, [], .7);
    title(i)
    shg
 end
