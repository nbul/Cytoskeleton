%Folders with images
cent_dir_control =[filedir, '/control/centrosomes'];
b_dir_control = [filedir, '/control/borders'];
cent_dir_mutant =[filedir, '/mutant/centrosomes'];
b_dir_mutant = [filedir, '/mutant/borders'];

%Folder to save information about cells
if exist([filedir,'/control/images_analysed'],'dir') == 0
    mkdir(filedir,'/control/images_analysed');
end
im_dir_control = [filedir, '/control/images_analysed'];

if exist([filedir,'/control/summary'],'dir') == 0
    mkdir(filedir,'/control/summary');
end
sum_dir_control = [filedir, '/control/summary'];

%Folder to save information about cells for mutant
if exist([filedir,'/mutant/images_analysed'],'dir') == 0
    mkdir(filedir,'/mutant/images_analysed');
end
im_dir_mutant = [filedir, '/mutant/images_analysed'];

if exist([filedir,'/mutant/summary'],'dir') == 0
    mkdir(filedir,'/mutant/summary');
end
sum_dir_mutant = [filedir, '/mutant/summary'];

if exist([filedir,'/mutant/summary/cellbycell'],'dir') == 0
    mkdir(filedir,'/mutant/summary/cellbycell');
end
sum_dir_mutant_cs = [filedir, '/mutant/summary/cellbycell'];

if exist([filedir,'/control/summary/cellbycell'],'dir') == 0
    mkdir(filedir,'/control/summary/cellbycell');
end
sum_dir_control_cs = [filedir, '/control/summary/cellbycell'];