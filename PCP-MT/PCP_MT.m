% Script for Automatic Analysis of N embryos
% The number of filesets should be entered for each experiment
% Adapted to use Sobel 5x5 to find direction
% Check analysis of Microtubule density

clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
dens_dir =[filedir, '/cytoskeleton_average'];
MTSD_dir =[filedir, '/cytoskeleton'];
b_dir = [filedir, '/borders'];

usedefault = questdlg(strcat('Is there PCA data'),'Settings','Yes','No','No');
if strcmp(usedefault, 'No')
    choice = 1;
else
    choice = 0;
end

%Folder to save information about cells
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];

%Folder to save information about cells
if exist([filedir,'/images_MTs'],'dir') == 0
    mkdir(filedir,'/images_MTs');
end
im_dir2 = [filedir, '/images_MTs'];

%Folder to save information about cells
if exist([filedir,'/images_Cells'],'dir') == 0
    mkdir(filedir,'/images_Cells');
end
im_dir3 = [filedir, '/images_Cells'];

if choice == 0
    %Folder to save information about cells
    if exist([filedir,'/images_PCP'],'dir') == 0
        mkdir(filedir,'/images_PCP');
    end
    im_dir4 = [filedir, '/images_PCP'];
    pcp_dir = [filedir, '/pcp_polarity'];
end

if exist([filedir,'/summary'],'dir') == 0
    mkdir(filedir,'/summary');
end
sum_dir = [filedir, '/summary'];

if exist([filedir, '/distribution'],'dir') == 0
    mkdir(filedir,'/distribution');
end
dist_dir = [filedir, '/distribution'];



%% Number of files to analyse
cd(dens_dir);
files = dir('*.tif');
cd(currdir);

Averages = zeros(length(files),7);

if choice == 0
    fractions = zeros(length(files), 12);
else
    fractions = zeros(length(files), 8);
end

allmu = 0;
allor = 0;
allPCP = 0;
allmu2 = 0;
allor2 = 0;
allPCP2 = 0;
headers = {'Embryo', 'Signal area', 'sem','Signal intensity','sem',...
    'Area', 'sem','Eccentricity','sem','Direction_cell','sem', ...
    'SD', 'sem','DEV','sem', 'Elongation','sem', 'Alignment','sem',...
    'Cell number','Outliers'};

Averages_filename = 'Summary.csv';
if exist([filedir,'/summary/cellbycell'],'dir') == 0
    mkdir(filedir,'/summary/cellbycell');
end
result_dir = [filedir, '/summary/cellbycell'];

%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';

for loop=1:length(files)
    %% reading files
    clear Name Number Actin_file Image_actin Path  Image_borders
    cd(dens_dir);
    Actin_file = [num2str(loop),'.tif'];
    Image = imread(Actin_file);
    cd(MTSD_dir);
    Actin_file = [num2str(loop),'.tif'];
    Image2 = imread(Actin_file);
    cd(b_dir);
    Path = [b_dir, '/', num2str(loop),'/'];
    cd(Path);
    Image_borders = imread('handCorrection.tif');
    [im_x, im_y] = size(Image);
    
    if choice == 0
        Path = [pcp_dir, '/', num2str(loop),'/Result/'];
        cd(Path);
        PCPdata = csvread('PCA_Cell-by-Cell_Polarity.csv', 1, 0);
        PCPangle = PCPdata(:,3);
    end
    
    %% Collect data about cells and boundaries
    borders;
    
    %% Cell-by-cell analysis
    
    cellbycelldensity;    
    
    MTSD;
    vonmises;
    
    %% writing down summarised data
    summaries;
    
    %% Open MTs image, adjust and generate cell masks
    imageprocessing;
    
end

allmu(1) = [];
allor(1) = [];
allmu2(1) = [];
allor2(1) = [];

cd(sum_dir);
 %% replace in the next line allmu with allmu2 to plot non-normalised distribution
figure1 = polarhistogram(deg2rad([allmu; allmu + 180]), 18, 'FaceColor', '#77AC30', 'LineWidth',2);
set(gca,'FontSize',36);
rticklabels([])
saveas(figure1, 'MT_orientation.eps', 'epsc');

%% replace in the next line allor with allor2 to plot non-normalised distribution
figure2 = polarhistogram(deg2rad([allor; allor + 180]), 18, 'FaceColor', '#777777', 'LineWidth',2);
set(gca,'FontSize',36);
rticklabels([])
saveas(figure2, 'Cell_orientation.eps', 'epsc');

if choice == 0
    allPCP(1) = [];
    allPCP2(1) = [];
    %% replace in the next line allPCP with allPCP2 to plot non-normalised distribution
    figure2 = polarhistogram(deg2rad([allPCP; allPCP + 180]), 18, 'FaceColor', '#00FFFF', 'LineWidth',2);
    set(gca,'FontSize',36);
    rticklabels([])
    saveas(figure2, 'PCP_orientation.eps', 'epsc'); 
end

% writetable(Diffangle, 'all_angles.csv');
Averages = sortrows(Averages,1);

csvwrite_with_headers(Averages_filename,Averages,headers);


%% comparing medians
[pval3, med3, P3] = circ_cmtest(deg2rad(allor2),deg2rad(allmu2));
[pval4, med4, P4] = circ_cmtest(deg2rad(allor),deg2rad(allmu));

if choice == 0
    [pval1, med1, P1] = circ_cmtest(deg2rad(allor2),deg2rad(allPCP2));
    [pval2, med2, P2] = circ_cmtest(deg2rad(allor),deg2rad(allPCP));
    [pval5, med5, P5] = circ_cmtest(deg2rad(allPCP2),deg2rad(allmu2));
    [pval6, med6, P6] = circ_cmtest(deg2rad(allPCP),deg2rad(allmu));
    
    headers2 = {'cell-PCP', 'cell-PCP-norm','cell-MT','cell-MT-norm','PCP-MT','PCP-MT-norm'};
    csvwrite_with_headers('pvalues.csv',[pval1,pval2,pval3, pval4, pval5, pval6],headers2);
    headers3 = {'CellPD','CellAP','CellPD_norm','CellAP_norm','MTPD','MTAP',...
        'MTPDnorm','MTAPnorm','PCPPD','PCPAP', 'PCPPDnorm','PCPAPnorm'};
    csvwrite_with_headers('fractions.csv',fractions,headers3);
else
    headers2 = {'cell-MT','cell-MT-norm'};
    csvwrite_with_headers('pvalues.csv',[pval3, pval4],headers2);
    headers3 = {'CellPD','CellAP','CellPD_norm','CellAP_norm','MTPD','MTAP','MTPDnorm','MTAPnorm'};
    csvwrite_with_headers('fractions.csv',fractions,headers3);
end

cd(currdir);

clc
clear variables
close all
