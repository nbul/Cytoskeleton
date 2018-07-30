%% Writing down summarised data for MTSD
cd(MTSD_res_dir);

summary_MTSD = [cell_data(:,1),cell_data(:,3:5),SD',mu',...
    1./sqrt(1-cell_data(:,4).*cell_data(:,4)), 100*(erf(10./SD'/sqrt(2))-erf(-10./SD'/sqrt(2)))/2];
summary_MTSD(summary_MTSD(:,6)<0,6) = summary_MTSD(summary_MTSD(:,6)<0,6) + 180;

summary_filename = [num2str(Number),'_MTSD.csv'];
csvwrite_with_headers(summary_filename,summary_MTSD,headers_MTSD_ind);



%% Obtaining averaged data for MTSD
Averages_MTSD(loop,1) = Number;
Averages_MTSD(loop,2:2:size(summary_MTSD(:,2:end),2)*2) = mean(summary_MTSD(:,2:end),1);
Averages_MTSD(loop,3:2:(size(summary_MTSD(:,2:end),2)*2+1))...
    = sqrt(var(summary_MTSD(:,2:end),1)/size(summary_MTSD(:,2:end),1));

%% Writing down summarised data for density
cd(dens_res_dir);

summary_dens = [cell_data(:,1),intensity',intensity_mts', mts_area', mts_density'...
    mts_bundling', Uniformity', UNAAD', Spars', skew', kurt', space',...
    Ent1',Ent2', Sdr', Sdq', cell_data(:,3:4)];

outlier = isoutlier(summary_dens(:,end-1))+isoutlier(summary_dens(:,end));
summary_dens(outlier > 0,:) = [];

summary_filename = [num2str(Number),'_density.csv'];
csvwrite_with_headers(summary_filename,summary_dens,headers_dens_ind);

%% Obtaining averaged data for density
Averages_dens(loop,1) = Number;
Averages_dens(loop,2:2:size(summary_dens(:,2:end),2)*2) = mean(summary_dens(:,2:end),1);
Averages_dens(loop,3:2:(size(summary_dens(:,2:end),2)*2+1))...
    = sqrt(var(summary_dens(:,2:end),1)/size(summary_dens(:,2:end),1));
Averages_dens(loop,size(summary_dens(:,2:end),2)*2+2) = size(summary_dens(:,2:end),1);
Averages_dens(loop,size(summary_dens(:,2:end),2)*2+3) = sum(outlier > 0);
cd(currdir);