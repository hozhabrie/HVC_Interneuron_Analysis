function vigiSpecGUI(varargin)
%input (in that order!)
%   ax: axes to plot on
%   w: signal (1 dimensional vector)
%   fs: sample freq (hz)
%   F: vector of frequencies to sample at (see defaultArgs)
%   thr: 0-1 scaling parameter to raise the noise threshold
%   window: 2^x points for fft (see defaultArgs)
%   step: sliding window parameter (= overlap of samples between window n and n+1)
%   tRange: if you want to shift your time vector
%   freeze colors: this is an advanced option if you want to plot something
%   on top of the spectrogram without affecting the spectrogram itself
defaultArgs= {gca,[], [], 500:5:8e3,.6,512,384,[],0};%this is standard for bird
optargs=defaultArgs;
optargs(1:nargin)=varargin(1:nargin);
optargs(cellfun(@isempty,optargs))=defaultArgs(cellfun(@isempty,optargs));
[ax,w,fs,F,thr,win,step,tR,rgb] = optargs{:};

    
[~,F,T,P]=spectrogram(w,win,step,F,fs);%make it
if ~isempty(tR)
    T=linspace(tR(1),tR(2),length(T));
end
S=10*log10(P);
S(S==-Inf)=min(isreal(S(:)));
% cmap=jet(256);    
% cmap(1:8,1:3)=zeros(8,3);    
% cmap=bone(256);
cmap=parula(256);    
cmap(1:8,1:3)=zeros(8,3);    
colormap(cmap);
if ~rgb
    try
        imagesc(ax,T,F,S);
        set(ax,'clim',[min(S(:))+thr*range(S(:)),max(S(:))]);%change the colors
        set(ax,'ydir','normal')
    catch
        axis(ax)
        imagesc(T,F,S);
        set(ax,'clim',[min(S(:))+thr*range(S(:)),max(S(:))]);%change the colors
        set(ax,'ydir','normal')
    end
else
    map=[min(S(:))+thr*range(S(:)),1.1*max(S(:))];
%     map=[min(S(:))+thr*range(S(:)),1.1*max(S(:))];
    S(S<map(1))=map(1);
    S(S>map(2))=map(2);
    S=ceil((S-map(1))/range(map)*255)+1;
    Srgb=ind2rgb(S,cmap);
    imagesc(ax,T,F,Srgb);set(gca,'ydir','normal')
end