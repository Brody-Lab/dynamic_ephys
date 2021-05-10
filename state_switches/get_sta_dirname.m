function stadir = get_sta_dirname(p)
dp = set_dyn_path;
stalab = p.which_switch;
if sum(p.model_smooth_wdw)>1
    stalab = sprintf('%s_mvnwdw_%i',stalab,p.model_smooth_wdw);
end
if sum(p.t_buffers)>0
    stalab = sprintf('%s_tbuffs_%i_%i',stalab,p.t_buffers(1)*1e3,p.t_buffers(2)*1e3);
end

if sum(p.min_pre_dur + p.min_pre_dur) > 0
    stalab = sprintf('%s_pre%i_post%i', ...
        stalab,p.min_pre_dur*1e3,p.min_post_dur*1e3);
end
if p.clear_bad_strengths
    stalab = sprintf('%s_bad%i', stalab, p.bad_strength);
end
stadir = fullfile(dp.sta_dir,stalab);
if ~exist(stadir,'dir'), mkdir(stadir); end