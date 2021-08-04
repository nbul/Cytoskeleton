%% A script to calculate overall alignment between the cell major axis and main direction of the MT network

clc;
close all;
clear variables;

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

celldir = [filedir, '/summary/cellbycell'];
sumdir = [filedir, '/summary'];

% answer = inputdlg({'Enter eccentricity','Enter error'},'Input',1,{'0.75','0.025'});
% ecc = str2double(answer(1));
% error = str2double(answer(2));
cd(celldir);
files = dir('*.csv');

%% reading data

celldata = zeros(1,6);
for loop = 1:numel(files)
    cd(celldir);
    celldata = [celldata; csvread(['Summary_MTSD',num2str(loop),'.csv'], 1,4)]; %#ok<AGROW>
end

celldata(1,:) = [];

average= zeros(4,6);
counter2 = 0;
error = 0.025;
for ecc=0.7:0.05:0.85
    %% Leaving only cells of the correct eccentricity
    counter = 0;
    counter2 = counter2 + 1;
    celldataclean = zeros(1,4);
    for i = 1:size(celldata,1)
        if celldata(i,1)>= (ecc - error) && celldata(i,1) < (ecc + error)
            counter = counter + 1;
            celldataclean(counter,1) = celldata(i,1);
            celldataclean(counter,2) = celldata(i,4);
            celldataclean(counter,3) = celldata(i,2);
            celldataclean(counter,4) = abs(celldata(i,4)-celldata(i,2));
            if celldataclean(counter,4) > 90
                celldataclean(counter,4) = 180 - celldataclean(counter,4);
            end
        end
    end
    
    [dataout, TF] = rmoutliers(celldataclean(:,4), 'quartiles');
    average(counter2,:) = [mean(celldataclean(:,1),1), std(celldataclean(:,1),0,1), mean(dataout,1), std(dataout), counter, sum(TF)];
end

all = array2table(average);
all.Properties.VariableNames = {'Eccentricity', 'EccSD', 'CellMTalignment','CellMTSD', 'Ncells', 'Noutliers'};
cd(sumdir);
writetable(all,'Cell_MT_alignment.csv');
cd(currdir);

close all;
clear variables;