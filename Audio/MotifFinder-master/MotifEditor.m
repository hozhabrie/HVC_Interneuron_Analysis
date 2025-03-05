function varargout = MotifEditor(varargin)
% MOTIFEDITOR MATLAB code for MotifEditor.fig
%      MOTIFEDITOR, by itself, creates a new MOTIFEDITOR or raises the
%      existing==
%      singleton*.
%
%      H = MOTIFEDITOR returns the handle to a new MOTIFEDITOR or the handle to
%      the existing singleton*.
%
%      MOTIFEDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTIFEDITOR.M with the given input arguments.
%
%      MOTIFEDITOR('Property','Value',...) creates a new MOTIFEDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MotifEditor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MotifEditor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MotifEditor

% Last Modified by GUIDE v2.5 11-Oct-2017 12:19:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MotifEditor_OpeningFcn, ...
                   'gui_OutputFcn',  @MotifEditor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MotifEditor is made visible.
function MotifEditor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MotifEditor (see VARARGIN)

% Choose default command line output for MotifEditor
handles.output = hObject;
% addpath fcns%add it to path beforehand
linkaxes([handles.axes1,handles.axes2],'x');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MotifEditor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MotifEditor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pullupdata(hObject);
plotdata(hObject,1);
axis tight

% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
% hObject    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Motif=handles.Motif;
f=handles.currF;
plotdata(hObject,0);
x=ginput(2);
loc=x(:,1);
loc = sort(loc);
Motif(f).start=[Motif(f).start,loc(1)];
Motif(f).stop=[Motif(f).stop,loc(2)];
Motif(f).center=[Motif(f).center,mean(loc)];
Motif(f).warp=[Motif(f).warp,range(loc)/handles.tempLength];
[~,inds]=sort(Motif(f).center);
Motif(f).start=Motif(f).start(inds);
Motif(f).stop=Motif(f).stop(inds);
Motif(f).center=Motif(f).center(inds);
Motif(f).warp=Motif(f).warp(inds);
handles.Motif=Motif;
guidata(hObject,handles)
plotdata(hObject,0);

% --- Executes on button press in move.
function move_Callback(hObject, eventdata, handles)
% hObject    handle to move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Motif=handles.Motif;
f=handles.currF;
x=ginput(1);
loc=x(1);
[startD,indStart]=min(abs(Motif(f).start-loc));
[stopD,indStop]=min(abs(Motif(f).stop-loc));

if startD<stopD
    Motif(f).start(indStart)=NaN;
    handles.Motif=Motif;
    guidata(hObject,handles)
    plotdata(hObject,0);
    x=ginput(1);
    Motif(f).start(indStart)=x(1);
    Motif(f).warp(indStart)=(Motif(f).start(indStart)-Motif(f).stop(indStart))/handles.tempLength;
    Motif(f).center(indStart)=mean([Motif(f).start(indStart),Motif(f).stop(indStart)]);
end
if stopD<startD
    Motif(f).stop(indStop)=NaN;
    handles.Motif=Motif;
    guidata(hObject,handles)
    plotdata(hObject,0);
    x=ginput(1);
    Motif(f).stop(indStop)=x(1);
    Motif(f).warp(indStop)=(Motif(f).start(indStop)-Motif(f).stop(indStop))/handles.tempLength;
    Motif(f).center(indStop)=mean([Motif(f).start(indStop),Motif(f).stop(indStop)]);
end
handles.Motif=Motif;
guidata(hObject,handles)
plotdata(hObject,0);

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Motif=handles.Motif;
f=handles.currF;
x=ginput(1);
loc=x(1);
[~,indRm]=min(abs(Motif(f).center-loc));
Motif(f).start(indRm)=[];
Motif(f).stop(indRm)=[];
Motif(f).center(indRm)=[];
Motif(f).warp(indRm)=[];
handles.Motif=Motif;
guidata(hObject,handles)
plotdata(hObject,0);

% --- Executes on button press in SummaryPlot.
function SummaryPlot_Callback(hObject, eventdata, handles)
% hObject    handle to SummaryPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Motif=handles.Motif;
figure(1);clf;hold on;
mi=0;
offset=4e4*.2;
loc=nan(1,length(Motif));
for t=1:length(Motif)
    wave=audioread([handles.F,handles.files{t}]);
    ms=round(Motif(t).start*handles.fs);
    me=round(Motif(t).stop*handles.fs);
    o=round(Motif(t).center*handles.fs);
    w=Motif(t).warp;
    startM=mi;
    for m=1:length(ms)
        start=ms(m)-offset;
        stop=me(m)+offset;
        audio=wave(start:stop);
        tN=(1:length(audio))+ms(m)-o(m)-offset;%center it
        tN=tN/4e4/w(m);
        if t==handles.currF
            plot(tN,audio+mi,'r')
        else
            plot(tN,audio+mi,'k')
        end
        line(tN(offset)*[1,1],mi+[-1,1])
        line(tN(end-offset)*[1,1],mi+[-1,1])
        mi=mi+1;
    end
    loc(t)=(startM+mi-1)/2;
    mi=mi+2;
    line((handles.tempLength/2+offset/4e4)*[-1,1],(mi-1)*[1,1],'linestyle',':')
end
[wave, fs] = audioread(handles.template);
shift = length(wave)/fs/2;
plot((1:length(wave))/fs - shift, wave/max(abs(wave))*3 + mi + 1, 'b');
axis tight
set(gca,'ytick',[loc mi + 1]);
set(gca,'yticklabel',[num2strVec(1:t); 'Template']);
ylabel('Trial number');
xlabel('Time (s)');
% --- Executes on button press in chngThr.
function chngThr_Callback(hObject, eventdata, handles)
% hObject    handle to chngThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Motif = handles.Motif;
prompt = {'Threshold'};
dlg_title = 'Change Threshold';
% defaultans = {num2str(Motif(handles.currF).thres)};
defaultans = {num2str(0.9)};
answer = inputdlg(prompt,dlg_title,[1,50],defaultans);
thresh=str2num(answer{1});
update = 0;
if thresh>0&&thresh<1
    [start,stop,center,warp]=findMotifs([handles.F,handles.files{handles.currF}],handles.template,thresh);
    numberOfMotifsFound = length(start);
    if numberOfMotifsFound
        Motif(handles.currF).start = start;
        Motif(handles.currF).stop = stop;
        Motif(handles.currF).center = center;
        Motif(handles.currF).warp = warp;
        Motif(handles.currF).thresh = thresh;
        handles.Motif = Motif;
        update = 1;
    else
        msgbox('Please try again with lower threshold... if there are motifs here...')
    end
else
    error('Last time I checked, correlations needed to be between 0 and 1')
end
guidata(hObject,handles)
if update
    plotdata(hObject,0);
end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
Motif=handles.Motif;
save([handles.F,'MotifTimes.mat'],'Motif')
msgbox('Jesus saves, and so do you!!')

% --- Executes on button press in prev.
function prev_Callback(hObject, eventdata, handles)
% hObject    handle to prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.currF=max(1,handles.currF-1);
guidata(hObject,handles);
plotdata(hObject,1);
axis(handles.axes1,'tight');


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
handles.currF=min(length(handles.files),handles.currF+1);
guidata(hObject,handles);
plotdata(hObject,1);
axis(handles.axes1,'tight');

% --- Executes on slider movement.
function timeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to timeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles=guidata(hObject);
center=get(hObject,'Value') /(get(hObject,'Max') -get(hObject,'Min') );
currR=range(xlim(handles.axes1));
info=audioinfo([handles.F,handles.files{handles.currF}]);
tRange=info.Duration-currR;
newX=center*tRange+currR/2;
xlim(handles.axes1,newX+[-1,1]*currR/2);

% --- Executes during object creation, after setting all properties.
function timeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles=guidata(hObject);
f=str2double(get(hObject,'String'));
if rem(f,1)==0&&f>0&&f<=length(handles.files)
    handles.currF=f;
    guidata(hObject);
    guidata(hObject,handles);
    plotdata(hObject,1);
    axis(handles.axes1,'tight');
else
    msgbox('Have you tried giving a normal number maybe?')
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddStart.
function AddStart_Callback(hObject, eventdata, handles)
% hObject    handle to AddStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
Motif=handles.Motif;
f=handles.currF;
plotdata(hObject,0);
x=ginput(1);
loc=x(:,1);
loc = sort(loc);
Motif(f).start=[Motif(f).start,loc(1)];
Motif(f).stop=[Motif(f).stop,loc(1)+handles.tempLength];
Motif(f).center=[Motif(f).center,loc(1)+handles.tempLength/2];
Motif(f).warp=[Motif(f).warp,1];
[~,inds]=sort(Motif(f).center);
Motif(f).start=Motif(f).start(inds);
Motif(f).stop=Motif(f).stop(inds);
Motif(f).center=Motif(f).center(inds);
Motif(f).warp=Motif(f).warp(inds);
handles.Motif=Motif;
guidata(hObject,handles)
plotdata(hObject,0);