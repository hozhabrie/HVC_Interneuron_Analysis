function plotdata(hObject,specToo)
handles=guidata(hObject);
Motif=handles.Motif;
f=handles.currF;
[wave,fs]=audioread([handles.F,handles.files{f}]);
if specToo
    f=handles.currF;
    vigiSpecGUI(handles.axes1,wave,fs,500:5:8e3,.6);
    set(handles.TrialN,'String',[handles.files{handles.currF} '(' num2str(f) '/' num2str(length(handles.files)) ')'])
else
    xr=xlim;
end

prevLines= findobj(handles.axes1, 'Type', 'line');
if ~isempty(prevLines)
    delete(prevLines);
end
prevLines= findobj(handles.axes2, 'Type', 'line');
if ~isempty(prevLines)
    delete(prevLines);
end

plot(handles.axes2,(1:length(wave))/fs,wave, 'k')
ylim(handles.axes2,[-1,1])



for m=1:length(Motif(f).start)
    line(handles.axes1,Motif(f).start(m)*[1,1],[500,8e3],'color','g', 'linewidth', 2)
    line(handles.axes1,Motif(f).stop(m)*[1,1],[500,8e3],'color','r', 'linewidth', 2)    
    line(handles.axes1,Motif(f).center(m)*[1,1],[500,8e3],'color','c', 'linewidth', 2)

    line(handles.axes2,Motif(f).start(m)*[1,1],[-1,1],'color','g', 'linewidth', 2)
    line(handles.axes2,Motif(f).stop(m)*[1,1],[-1,1],'color','r', 'linewidth', 2)    
    line(handles.axes2,Motif(f).center(m)*[1,1],[-1,1],'color','c', 'linewidth', 2)
end
if ~specToo
    xlim(handles.axes2,xr);
end