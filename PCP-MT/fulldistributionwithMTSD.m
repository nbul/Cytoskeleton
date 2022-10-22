clc;
close all;
clear variables;

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

distdir = [filedir, '/distribution'];
celldir = [filedir, '/summary/cellbycell'];
sumdir = [filedir, '/summary'];

cd(celldir);
files = dir('*.csv');

%% reading data

celldata = zeros(1,6);
distdata = zeros(45,1);
for loop = 1:numel(files)
    cd(celldir);
    celldata = [celldata; csvread(['Summary_MTSD',num2str(loop),'.csv'], 1,4)]; %#ok<AGROW>
    cd(distdir);
    distdata = [distdata, csvread([num2str(loop),'_distribution.csv'],0,1)]; %#ok<AGROW>
end

celldata(1,:) = [];
distdata(:,1) = csvread([num2str(loop),'_distribution.csv'],0,0,[0,0,44,0]);

MTSDdata = zeros(1,5);
error = 0.025;
for ecc=0.65:0.05:0.85
    
    %% Leaving only cells of the correct eccentricity
    counter = 0;
    celldataclean = zeros(1,3);
    distdataclean = distdata(:,1);
    for i = 1:size(celldata,1)
        if celldata(i,1)>= (ecc - error) && celldata(i,1) < (ecc + error)
            counter = counter + 1;
            celldataclean(counter,1) = celldata(i,1);
            celldataclean(counter,2) = celldata(i,4);
            celldataclean(counter,3) = celldata(i,2);
            distdataclean(:,counter+1) = distdata(:,i+1);
            
        end
    end
    
    DEV = celldataclean(:,2) - celldataclean(:,3);
    %% centering distributions
    distdatashifted = distdata(:,1);
    for i = 1:counter
        if celldataclean(i,2) >= 90
            celldataclean(i,2) = celldataclean(i,2) - 180;
        end
        counter2 = 0;
        for k = 1:45
            if  distdataclean(k,1) >= celldataclean(i,2)
                counter2 = counter2 + 1;
                distdatashifted(counter2,i+1) = distdataclean(k,i+1);
            end
        end
        for k = 1:45
            if  distdataclean(k,1) < celldataclean(i,2)
                counter2 = counter2 + 1;
                distdatashifted(counter2,i+1) = distdataclean(k,i+1);
            end
        end
    end
    %% shifting distributions
    for i = 1:45
        if distdatashifted(i,1) >= 0
            distdatashifted(i,1) = distdatashifted(i,1) - 88;
        else
            distdatashifted(i,1) = distdatashifted(i,1) + 92;
        end
    end
    
    distdatashifted = sortrows(distdatashifted,1);
    
    %% average statistics
    finaldata = zeros(1,4);
    
    finaldata(1,1) = mean(celldataclean(:,1),1); % mean eccentricity
    finaldata(1,2) = std(celldataclean(:,1),0,1); % stardanrd deviation
    finaldata(1,3) = finaldata(1,2)/sqrt(counter - 1); % standard error
    finaldata(1,4) = counter; % number of cells
    
    finaldata(3:47,1) = distdata(:,1); % angle bins
    finaldata(3:47,2) = mean(distdatashifted(:,2:end),2); % average distribution
    finaldata(3:47,3) = std(distdatashifted(:,2:end),0,2); % standard deviation
    finaldata(3:47,4) = finaldata(3:47,3)/sqrt(counter - 1); % standard error
    
    edges2 = ceil(min(DEV)):1:ceil(max(DEV));
    [N,edges] = histcounts(DEV,edges2);
    bincenters = edges(1:end-1)+0.5;
    %% writing file
    cd(sumdir);
    csvwrite(['ecc_', num2str(ecc),'_error_',num2str(error),'.csv'], finaldata);
    csvwrite(['ecc_', num2str(ecc),'_error_',num2str(error),'_DEV.csv'], [bincenters', N']);
    
    m_added_norm = distdataclean;
    vonmises;
    
    %datatemp = rmoutliers([celldataclean(:,1),SD'],1);
    datatemp = [celldataclean(:,1),SD'];
    
    MTSDdata = [MTSDdata; mean(datatemp(:,1)), std(datatemp(:,1)), ...
        mean(datatemp(:,2)), std(datatemp(:,2)), length(datatemp(:,2))];
    
end

MTSDdata(1,:) = [];
headers = {'ECC', 'ECC_SD', 'MTSD', 'MTSD_SD', 'Cells'};
csvwrite_with_headers('Ecc_vs_MTSD.csv',MTSDdata,headers);
cd(currdir);

close all;
clear variables;