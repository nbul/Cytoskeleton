clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

files = dir('*.csv');

number = numel(files);

data_all = zeros(1,6);
stats = zeros(1,31);
for i=1:number
    
    data = readtable([num2str(i),' Spots.csv']);
    tracks = unique(data.TRACK_ID);
    deltaX = zeros(1, length(tracks));
    deltaY = zeros(1, length(tracks));
    angle = zeros(1, length(tracks));
    speed = zeros(1, length(tracks));
    duration = zeros(1, length(tracks));
    distance = zeros(1, length(tracks));
    
    for k = 1: length(tracks)
        T = data(data.TRACK_ID == tracks(k),:);
        deltaX(k) = T.POSITION_X(end) - T.POSITION_X(1);
        deltaY(k) = -(T.POSITION_Y(end) - T.POSITION_Y(1));
        angle(k) = atan2(deltaY(k), deltaX(k)) * 180 / pi;
        duration(k) = (size(T,1) - 1) * 1.116; % in seconds
        speed(k) = sqrt(deltaX(k) * deltaX(k) + deltaY(k) * deltaY(k)) * 0.09 /duration(k); % in microns
        distance(k) = sqrt(deltaX(k) * deltaX(k) + deltaY(k) * deltaY(k)) * 0.09;
    end
    
    data_movie = array2table([tracks, angle', duration', speed', distance']);
    data_movie.Properties.VariableNames = {'TRACK_ID', 'ANGLE', 'DURATION', 'SPEED', 'DISTANCE'};
    writetable(data_movie,[num2str(i),' Track result.csv']);
    
    data_all = [data_all; [ones(length(tracks),1) * i, tracks, angle', duration', speed', distance']]; %#ok<AGROW>
    
    
    stats(i,1) = i;
    stats(i,2) = length(tracks);
    
    stats(i,3) = mean(data_movie.DURATION);
    stats(i,4) = std(data_movie.DURATION);
    
    stats(i,5) = mean(data_movie.SPEED);
    stats(i,6) = std(data_movie.SPEED);
    
    stats(i,7) = mean(data_movie.DISTANCE);
    stats(i,8) = std(data_movie.DISTANCE);
    
    % 0-30°
    stats(i,9) = size(data_movie(abs(sind(data_movie.ANGLE))<0.5,:),1)/length(tracks)*100; % percent tracks with angles less then 30 with X-axis;
    stats(i,10) = mean(data_movie(abs(sind(data_movie.ANGLE))<0.5,:).DURATION); % duration for angles less then 30 with X-axis;
    stats(i,11) = std(data_movie(abs(sind(data_movie.ANGLE))<0.5,:).DURATION);
    
    stats(i,12) = mean(data_movie(abs(sind(data_movie.ANGLE))<0.5,:).SPEED); % speed for angles less then 30 with X-axis;
    stats(i,13) = std(data_movie(abs(sind(data_movie.ANGLE))<0.5,:).SPEED);
    
    stats(i,14) = mean(data_movie(abs(sind(data_movie.ANGLE))<0.5,:).DISTANCE); % distance for angles less then 30 with X-axis;
    stats(i,15) = std(data_movie(abs(sind(data_movie.ANGLE))<0.5,:).DISTANCE);
    
    % 30-60%
    stats(i,16) = size(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:),1)/length(tracks)*100; % percent tracks with angles greater then 30 but less then 60 with X-axis;
    stats(i,17) = mean(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:).DURATION); % duration for angles greater then 30 but less then 60 with X-axis;
    stats(i,18) = std(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:).DURATION);
    
    stats(i,19) = mean(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:).SPEED); % speed for angles greater then 30 but less then 60 with X-axis;
    stats(i,20) = std(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:).SPEED);
    
    stats(i,21) = mean(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:).DISTANCE); % distance for angles greater then 30 but less then 60 with X-axis;
    stats(i,22) = std(data_movie(abs(sind(data_movie.ANGLE))>0.5 & abs(sind(data_movie.ANGLE))<sind(60),:).DISTANCE);
    
    %60-90
    stats(i,23) = size(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:),1)/length(tracks)*100; % percent tracks with angles greater than 60 with X-axis;
    stats(i,24) = mean(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:).DURATION); % duration for angles greater than 60 with X-axis;
    stats(i,25) = std(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:).DURATION);
    
    stats(i,26) = mean(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:).SPEED); % speed for angles greater than 60 with X-axis;
    stats(i,27) = std(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:).SPEED);
    
    stats(i,28) = mean(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:).DISTANCE); % distance for angles greater than 60 with X-axis;
    stats(i,29) = std(data_movie(abs(sind(data_movie.ANGLE))>sind(60),:).DISTANCE);
    
    %directionality
    stats(i,30) = size(data_movie(data_movie.ANGLE>-90 & data_movie.ANGLE<=90,:),1)/length(tracks)*100;  % posterior direction
    stats(i,31) = size(data_movie(data_movie.ANGLE>0 & data_movie.ANGLE<=180,:),1)/length(tracks)*100;  % up
    
end

data_all(1,:) = [];
data_all = array2table(data_all);
data_all.Properties.VariableNames = {'MOVIE','TRACK_ID', 'ANGLE', 'DURATION', 'SPEED', 'DISTANCE'};
writetable(data_movie,'All tracks result.csv');

stats = array2table(stats);
stats.Properties.VariableNames = {'MOVIE', 'NTRACKS','DURATION', 'DURSTD', 'SPEED', 'SPSTD','DISTANCE', 'DISSTD',...
    'N0-30', 'DUR0-30', 'DUR0-30STD', 'SP0-30', 'SP0-30STD','DISTANCE0-30', 'DISR0-30STD',...
    'N30-60', 'DUR30-60', 'DUR30-60STD', 'SP30-60', 'SP30-60STD','DISTANCE30-60', 'DIS30-60STD',...
    'N60-90', 'DUR60-90', 'DUR60-90STD', 'SP60-90', 'SP60-90STD','DISTANCE60-90', 'DIS60-90STD',...
    'POSTERIOR', 'UP'};
writetable(stats,'Summary.csv');