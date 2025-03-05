function [EventStruct] = getEventTimes(audioPath, fs, eventTemplate, audioLabels, label)

% Description: 1. Uses annotated labels from Audacity (formatted to a table in Matlab);
%                   to save times for each type of event label that was tracked;
%                   annoted audio labels: 'sb' = song bout; 
%                                                  'sm' = song motif;  
%                                                  'pb' = pb bout;
%                                                  'pm' = isolate pb motif;
%                                                  'pmb' = pb motif within pb bout;
%                                                  'contaminated....' = stuff that can't be considered as spont... 
%                                                                              (but also can't be assessed for modulation)

% create structures for each type of event with starts, stops, centers, and warps
EventStruct.file = [audioPath 'audioFull.wav'];
EventStruct.start = audioLabels.start(audioLabels.label == label);
EventStruct.stop = audioLabels.stop(audioLabels.label == label);
eventLength = EventStruct.stop-EventStruct.start;
EventStruct.center = EventStruct.start+eventLength./2;
if ~strcmp(label, 'sb') 
    templateLength = length(eventTemplate)/fs;
    EventStruct.warp = eventLength./templateLength;
end
