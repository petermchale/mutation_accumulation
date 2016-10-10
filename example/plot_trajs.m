close all
clear all
format long 

% parameters
save_png = false;
annotate = true;
filename_extension = '_transitions';
trajNum = 6;
trajectory_type = 'successful';
png_file_name = 'traj.png';
linewidth = 3;
markersize = 10;
dir_name = 'large L/r0.51';  % 'large L/large_initial_tumor'; 

% load monte carlo trajectories
filename = [dir_name '/' trajectory_type '_trajectory' num2str(trajNum) '.dat' filename_extension];
trajectory_mc = load(filename); 
tt_mc = trajectory_mc(:,1); 
nn_mc = trajectory_mc(:,2:end);
nn_total_mc = sum(nn_mc,2);

% calculate deterministic trajectories
time_span = [0 tt_mc(end)];
nn_det_init = trajectory_mc(1,2:(end-1));
trajectory_det = calculate_det_traj(time_span,nn_det_init, [dir_name '/main.in']);

% plot trajectories
figure('color','white','renderer','zbuffer','units','inches','position',[1 1 8 7])
axes('fontsize',25,'linewidth',linewidth,'ticklength',[0.03 10]);
plot_trajectories_2(trajectory_mc, '-', linewidth, markersize);
plot_trajectories_2(trajectory_det, '--', linewidth, markersize);
set(gca,'ytick',[1e0 1e2 1e4 1e6 1e8 1e10 1e12]) 
xlim([0 tt_mc(end)]) 
set(gca,'tickdir','out')

% % plot total tumor size
% plot(tt_mc,nn_total_mc,'linestyle','s','linewidth',linewidth,'markersize',markersize)

% remove annotation, if requested
if ~annotate
    set(gca,'xtick',[],'ytick',[])
    legend off
end

% save as png file, if requested
if save_png 
	export_fig(png_file_name,'-m3')
end



