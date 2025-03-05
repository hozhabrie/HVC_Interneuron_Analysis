function makeWavs(folder,addT,threshold)
%This function will create the wav files for eGUI to select the syllables
%INPUT: list of ABF files and where they are located
%OUTPUT: it's main 
if ~exist([folder,'eguiWavs\'],'dir')
    mkdir([folder,'eguiWavs\'])
    mkdir([folder,'eguiWavs\template\'])
else
    files2rm=dir([folder,'eguiWavs\*.wav']);
    for f=1:length(files2rm)
        delete([folder,'eguiWavs\' files2rm{f}]);
    end
end
files=dir([folder,'ABF\*.abf']);%findall the files in the work space
files={files.name};
audioStart=zeros(1,length(files));
for i=1:length(files)
    [d,si,~]=abfload([folder,'ABF\',files{i}]); %loads abf the file
    glass=d(:,3);%3rd column is the visual acccess
    glass=glass/max(glass);
    lookat=glass>.2;%find all the times that the birds could see eachother
    stop=min((find(lookat,1,'last')+round(1/(si*1e-6))*addT),size(d,1));%the index of stoppping is addT seconds aftyer the last sample that they birds could see eachother
    lookat=find(lookat,1,'first'):stop;%indices of cutting, from the first time they can see eachother until stop
    wav=d(lookat,2);%only care about when the glass is open
    y=wav/max(abs(wav));%normalize wave file
    a=strtok(files{i},'.');%change from abf to wav
    name=[a,'.wav'];
    audiowrite([folder,'eguiWavs\',filesep,name],y,round(1/(si*1e-6)));%write it
    audioStart(i)=find(glass>.2,1,'first');%SAVE WJHERE YOU CUT IT
end
save([folder,'eguiWavs\Exp.mat'],'files','audioStart');
disp('done')

folder=[folder, 'eguiWavs\'];
threshold=num2str(threshold);
template=[folder, 'template\template.wav'];
save StandardPaths.mat folder template threshold