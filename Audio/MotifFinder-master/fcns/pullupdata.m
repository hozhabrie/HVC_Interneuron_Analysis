function pullupdata(hObject)
handles=guidata(hObject);
load StandardPaths.mat
prompt = {'Folder Name','Template File','Threshold'};
dlg_title = 'Load Data';
defaultans = {fWav,...                         %Change this to last path
              fTemplate,...
                threshold};
answer = inputdlg(prompt,dlg_title,[1,70],defaultans);
if ~isempty(answer)
    F=answer{1};
    fTemplate=answer{2};
    thresh=str2num(answer{3});
    if ~strcmp(F(end),'\')&&~strcmp(F(end),'/')
        F=[F,'\'];
    end
    files=dir([F,'*.wav']);
    if isempty(files)
        error('Have you tried giving me a folder with wav files in it YOU DUMKOPF?')
    end
    if ~strcmp(fTemplate(end-3:end),'.wav')
        error('maybe next time give me a wave file for a template, ya?')
    end
    files={files.name};

    fWav=F;
    threshold=answer{3};
    analyze=1;
    prevAnalysis=[F,'MotifTimes.mat'];
    if exist(prevAnalysis,'file')
        button = questdlg('You already analyzed this. Redo?','Ack!','Redo','Load Prev','Load Prev');
        if strcmp(button,'Load Prev')
            analyze=0;
        else
            analyze=1;
        end
    end

    if analyze
        for f=1:length(files)
            fprintf(['Analyzing file #' num2str(f) '/' num2str(length(files))])
            [start,stop,center,warp]=findMotifs([F,files{f}],fTemplate,thresh);
            Motif(f).file=files{f};
            Motif(f).start=start;
            Motif(f).stop=stop;
            Motif(f).center=center;
            Motif(f).warp=warp;
            Motif(f).thresh=thresh;
        end
        save([F,'MotifTimes.mat'],'Motif');
   %     save([F,'MotifTimes_FirstPass.mat'],'Motif');
    else
        load(prevAnalysis) 
    end
    
    info = audioinfo([F,files{1}]);
    handles.currF=1;
    handles.F=F;
    handles.fs=info.SampleRate;
    handles.files=files;
    handles.template=fTemplate;
    handles.Motif=Motif;
    info=audioinfo(fTemplate);
    handles.tempLength=info.Duration;
    save StandardPaths.mat threshold fWav fTemplate
    
    [wav,fs]=audioread(fTemplate);
    plot(handles.axes5,(1:length(wav))/fs,wav,'k');
    axis tight
    vigiSpecGUI(handles.axes4,wav,fs);
    
    guidata(hObject,handles);
else
    error('no file provided')
end