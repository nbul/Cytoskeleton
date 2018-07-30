dens_dir =[filedir, '/cytoskeleton_average'];
MTSD_dir =[filedir, '/cytoskeleton'];
b_dir = [filedir, '/borders'];
th_dir = [filedir, '/summary/threshold'];
sum_dir = [filedir, '/summary'];

%Folder to save information about cells
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];

if exist([filedir, '/distribution'],'dir') == 0
    mkdir(filedir,'/distribution');
end
dist_dir = [filedir, '/distribution'];

if exist([filedir,'/summary/density'],'dir') == 0
    mkdir(filedir,'/summary/density');
end
dens_res_dir = [filedir, '/summary/density'];

if exist([filedir,'/summary/MTSD'],'dir') == 0
    mkdir(filedir,'/summary/MTSD');
end
MTSD_res_dir = [filedir, '/summary/MTSD'];

headers_MTSD_ind = {'Cell', 'Area', 'Eccentricity','Direction_cell', ...
    'SD', 'DEV', 'Elongation', 'Alignment'};
headers_MTSD_all = {'Cell', 'Area', 'sem', 'Eccentricity', 'sem', 'Direction_cell',  'sem',...
    'SD', 'sem', 'DEV', 'sem', 'Elongation', 'sem', 'Alignment' 'sem',};

headers_dens_ind = {'Cell', 'Cell intensity', 'MTs intensity', 'Signal area', 'Density',...
    'Bundling', 'Uniformity', 'UNAAD', 'Sparseness', 'Skewness', 'Kurtosis',...
    'Gaps', 'Cell entropy', 'MTs entropy','Sdr','Sdq', 'Area', 'Eccentricity',};
headers_dens_all = {'Embryo', 'Cell intensity', 'sem', 'MTs intensity', 'sem',...
    'Signal area', 'sem','Density','sem','Bundling','sem',...
    'Uniformity','sem', 'UNAAD','sem','Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'Gaps', 'sem', 'Cell entropy', 'sem', 'MTs entropy', 'sem', 'Sdr', 'sem', 'Sdq', 'sem',...
    'Area', 'sem', 'Eccentricity','sem', 'Cell number', 'Outliers', };