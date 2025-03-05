function [starts,stops,centers,warps]=findMotifs(audioF,templateF,thresh)
% Input
%   audioF: string filename of .wav for a trial (~25 s)
%   templateF: string filename of a .wav file of the template song (~1 s)
%   thresh: threshold on the correlations,m between 0 and 1
% Output
%   starts: song onsets
%   stops, centers you know what to do 
%   warp: warp=length of song/length(template)

th = thresh; %threshold for accepting a motif (0-1)
ds=100; %downsampling factor
cutoutaddon=.2; % time in seconds added to template cutouts (at start and end) for warping
samplingrate=30000; %in Hz
cutoutaddon=cutoutaddon*samplingrate/ds;

motifTemplate_full = audioread(templateF); %modify template bounds if necessary, default: 1:end  %Vigi Data: (7819:26770)

%%%Run only to this point on the first run to cut out your template

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[audioMin]=audioBounds(motifTemplate_full,ds,1);
motifTemplate = audioMin;                       %final Motif template

%%%Creating a set of warped templates (default: 51 templates; 95 - 105 %)

k=1; %used for creating the warped templates
warping_range=[980:2:1020];
for f=warping_range    
        motifTemplate_warped=resample(motifTemplate,f,1000); %warpes template
        Warping.templates{k}=motifTemplate_warped;
        Warping.warp_factor{k}=f/1000;
    k=k+1;
end
clear k

MotifNum=0; %Index for writing MotifPositions

%%%Everything from this point needs to be done for every trial

CurrTrial_full=audioread(audioF);                %get current trial
[audioMin]=audioBounds(CurrTrial_full,ds,1);
CurrTrial=audioMin;


%%%Sliding Crosscorrelation
fprintf('\n find motifs')

offset = length(motifTemplate);
slidingCrossCorr = [];
ii=(length(CurrTrial)-offset);
for i = 1:ii
    tmp = corr(motifTemplate', CurrTrial(i:i+offset-1)');
    slidingCrossCorr = [slidingCrossCorr tmp];
    
    %waitbar(i / ii ) %%%slows loop down considerably 
end

%%%Finding peaks above threshold (th)

peakstart=[];
peakend=[];

for j=1:length(slidingCrossCorr)-1    
   
    if slidingCrossCorr(j)<=th && slidingCrossCorr(j+1)>=th
        peakstart=[peakstart j+offset/2];
    end
    
    if slidingCrossCorr(j)>=th && slidingCrossCorr(j+1)<=th
        peakend=[peakend j+offset/2];
    end
          
end
peaks=[peakstart' peakend'];


%%%now the peaks within threshold crossings and therefore the motifs are
%%%found
motifpos=[];

for y=1:size(peaks,1)
    [peak, nMax] = max(slidingCrossCorr(peaks(y,1)-offset/2:peaks(y,2)-offset/2));
    if y<3
    motifpos=[motifpos; peak nMax+peaks(y,1)-1];
    end
    
    if y>2  %this helps to avoid that several close crosscorr peaks within one motif cause many overlapping motif detections
        
        if nMax+peaks(y,1)-1 > motifpos(end,2)+cutoutaddon
        motifpos=[motifpos; peak nMax+peaks(y,1)-1];
        end
        
        if nMax+peaks(y,1)-1 < motifpos(end,2)+cutoutaddon
            
            if peak>motifpos(end,1)
                motifpos(end,:)=[];
                motifpos=[motifpos; peak nMax+peaks(y,1)-1];
            end
            
        end
    end 
  
end

numberOfMotifs = size(motifpos,1);
fprintf(['\n motifs detected: ' num2str(numberOfMotifs) '\n'])

if numberOfMotifs == 0
    fprintf(['Please try again with lower threshold... if there are motifs here...\n'])

    starts = [];
    stops = [];
    centers = [];
    warps = [];

else
    motifpos=motifpos(:,2)';


    motifonset=motifpos-offset/2;
    motifoffset=motifpos+offset/2;

    %%%now the best warping factor for each motif is found

    currmotifs=[];%needed for plotting (for indexing MotifPositions)

    %%%struct vectors for AllTrials
    warp_factors=[]; motif_centers=[]; motif_starts=[]; motif_ends=[];

    for e=1:size(motifpos,2)  %for each motif

        MotifNum=MotifNum+1;
        currmotifs=[currmotifs MotifNum];

        motifcutout=CurrTrial(motifpos(e)-offset/2-cutoutaddon:motifpos(e)+offset/2+cutoutaddon); %cuts out one motif +/- a little bit

        %figure  %this one is distributed (~13 lines and ~25 lines below)
        Crosscorr_maxpos=[];
        for f=1:size(Warping.templates,2) %for each warped template 
            warpoffset = length(Warping.templates{f});
            slidingCrossCorr_warp = [];
            for i = 1:(length(motifcutout)-warpoffset)  %sliding crosscorrelation
                tmp = corr(Warping.templates{f}', motifcutout(i:i+warpoffset-1)');
                slidingCrossCorr_warp = [slidingCrossCorr_warp tmp];
            end
            [maxcross, maxposcross]=max(slidingCrossCorr_warp);
            try
                %warping factor, crosscorr peak, position of crosscorr peak
                Crosscorr_maxpos=[Crosscorr_maxpos;...
                                  warping_range(f)/1000 maxcross maxposcross]; 
            catch
                keyboard
            end
    %                 ccAxis = (1:length(slidingCrossCorr_warp)) + warpoffset/2; %plotting the crosscorr for every warp factor
    %                 plot(ccAxis, slidingCrossCorr_warp);
    %                 hold on;
        end

        %save highest crosscorr peak and corresponding warp factor and position
        %within motifcutout

        MotifPositions(MotifNum).motifnum_local=e;

        MotifPositions(MotifNum).max_corr_val=max(Crosscorr_maxpos(:,2)); %max correlation value during warping
        [~,mwp]=max(Crosscorr_maxpos(:,2));
        MotifPositions(MotifNum).postwarp_shift=(Crosscorr_maxpos(mwp,3)+length(Warping.templates{mwp})/2)-ceil(length(motifcutout)/2); 
        MotifPositions(MotifNum).template_warp_factor=Crosscorr_maxpos(mwp,1);
        MotifPositions(MotifNum).postwarp_length=length(Warping.templates{mwp})*ds;
        MotifPositions(MotifNum).postwarp_pos=(motifpos(e)+MotifPositions(MotifNum).postwarp_shift)*ds;
        MotifPositions(MotifNum).postwarp_start=MotifPositions(MotifNum).postwarp_pos-MotifPositions(MotifNum).postwarp_length/2;
        MotifPositions(MotifNum).postwarp_end=MotifPositions(MotifNum).postwarp_pos+MotifPositions(MotifNum).postwarp_length/2;

        warp_factors=[warp_factors MotifPositions(MotifNum).template_warp_factor]; 
        motif_centers=[motif_centers MotifPositions(MotifNum).postwarp_pos]; 
        motif_starts=[motif_starts MotifPositions(MotifNum).postwarp_start]; 
        motif_ends=[motif_ends MotifPositions(MotifNum).postwarp_end];

    %     plot(motifcutout);
    %     line([ceil(length(motifcutout)/2), ceil(length(motifcutout)/2)],[0 1]); %plot prewarp center 
    %     line([ceil(length(motifcutout)/2)+MotifPositions(e).postwarp_shift,ceil(length(motifcutout)/2)+MotifPositions(e).postwarp_shift],[0 1],'color','k'); %plot postwarp center
    %     hold off     
    end

    starts = motif_starts/samplingrate;
    stops = motif_ends/samplingrate;
    centers = motif_centers/samplingrate;
    warps = warp_factors;
end

