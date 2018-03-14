%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(currdir);
javaaddpath(currdir);
filedir = uigetdir();
cd(filedir);

color0 = questdlg(strcat('Is there mCherry signal'),'Settings','Yes', 'No', 'Yes');
if strcmp(color0, 'Yes')
    CD8color = 2;
end

if strcmp(color0, 'Yes')
    color1 = questdlg(strcat('Which color is E-cad'),'Settings','Green', 'Blue','Green');
    color2 = questdlg(strcat('Which color is MTs'),'Settings','Green', 'Blue','Blue');
    
    if strcmp(color1, 'Blue')
        ColorCad = 0;
        ColorMT = 1;
    else
        ColorCad = 1;
        ColorMT = 0;
    end
        % folders for stripes
    if exist([filedir, '/CD8projection'],'dir') == 0
        mkdir(filedir, '/CD8projection');
    end
    stripe_dir = [filedir, '/CD8'];
else
    color1 = questdlg(strcat('Which color is E-cad'),'Settings','Green', 'Red', 'Blue','Green');
    color2 = questdlg(strcat('Which color is MTs'),'Settings','Green', 'Red', 'Blue','Blue');
    
    if strcmp(color1, 'Blue')
        ColorCad = 0;
        ColorMT = 1;
    elseif strcmp(color2, 'Blue')
        ColorCad = 1;
        ColorMT = 0;
    elseif strcmp(color1, 'Red') && strcmp(color2, 'Green')
        ColorCad = 1;
        ColorMT = 0;
    else
        ColorCad = 0;
        ColorMT = 1;
    end
end


%Folder with borders
if exist([filedir, '/borders'],'dir') == 0
    mkdir(filedir, '/borders');
end
border_dir = [filedir, '/borders'];

%Folder with cytoskeleton
if exist([filedir, '/cytoskeleton'],'dir') == 0
    mkdir(filedir, '/cytoskeleton');
end
cytM_dir = [filedir, '/cytoskeleton'];

%Folder with cytoskeleton averaged
if exist([filedir, '/cytoskeleton_average'],'dir') == 0
    mkdir(filedir, '/cytoskeleton_average');
end
cytA_dir = [filedir, '/cytoskeleton_average'];



% folders for originals
if exist([filedir, '/originals'],'dir') == 0
    mkdir(filedir, '/originals');
end
ori_dir = [filedir, '/originals'];


cd(filedir);
files = dir('*.czi');
%% Projections
for i=1:numel(files)
    name.original = [num2str(i), '_Out.czi'];
    original = bfopen(name.original);
    
    Series = original{1,1};
    seriesCount = size(Series, 1); %display size to check type of file
    Cyto.average = zeros(size(Series{1,1},1),size(Series{1,1},2));
    Cyto.max = zeros(size(Series{1,1},1),size(Series{1,1},2));
    CD8.average = zeros(size(Series{1,1},1),size(Series{1,1},2));
    Cad.max = zeros(size(Series{1,1},1),size(Series{1,1},2));
    
    if strcmp(color0, 'Yes')
        for plane = 1:(seriesCount/3)
            CD8.average = plus(CD8.average,double(Series{plane*3-CD8color,1}));
            Cad.max = max(Cad.max, double(Series{plane*3-ColorCad,1}));
            Cyto.max = max(Cyto.max, double(Series{plane*3-ColorMT,1}));
            Cyto.average = plus(Cyto.average,double(Series{plane*3-ColorMT,1}));
        end
        CD8.average = imadjust(uint16(CD8.average*3/seriesCount));
        Cyto.average = uint16(Cyto.average*3/seriesCount);
        Cad.background = imopen(Cad.max,strel('disk',50));
        Cad.max = Cad.max - Cad.background;
        Cad.max = imadjust(uint16(Cad.max));
        Cyto.max = uint16(Cyto.max);
        cd(stripe_dir);
        imwrite(CD8.average, [num2str(i), '.tif']);
    else
        for plane = 1:(seriesCount/2)
            CD8.average = plus(CD8.average,double(Series{plane*2-CD8color,1}));
            Cad.max = max(Cad.max, double(Series{plane*2-ColorCad,1}));
            Cyto.max = max(Cyto.max, double(Series{plane*2-ColorMT,1}));
            Cyto.average = plus(Cyto.average,double(Series{plane*2-ColorMT,1}));
        end
        Cyto.average = uint16(Cyto.average*2/seriesCount);
        Cad.background = imopen(Cad.max,strel('disk',50));
        Cad.max = Cad.max - Cad.background;
        Cad.max = imadjust(uint16(Cad.max));
        Cyto.max = uint16(Cyto.max);
    end
    
    cd(border_dir);
    imwrite(Cad.max, [num2str(i), '.tif']);
    cd(cytM_dir);
    imwrite(Cyto.max, [num2str(i), '.tif']);
    cd(cytA_dir);
    imwrite(Cyto.average, [num2str(i), '.tif']);
    cd(filedir);
    movefile([num2str(i), '_Out.czi'], ori_dir);
end

cd(currdir);
close all
clear variables
clc