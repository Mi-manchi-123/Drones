function out_csv = run_expert_multitarget_sweep()
%RUN_EXPERT_MULTITARGET_SWEEP Multi-target SNR and target-speed sweep.

this_dir = fileparts(mfilename('fullpath'));
if exist(fullfile(this_dir, 'ob_moco_v4_ablation_experiment.m'), 'file')
    scripts_dir = this_dir;
    root = this_dir;
    data_dir = fullfile(root, 'imu_replay');
    results_dir = fullfile(root, 'results_expert_multitarget_sweep');
else
    scripts_dir = this_dir;
    root = fileparts(this_dir);
    data_dir = fullfile(root, 'data', 'imu_replay');
    results_dir = fullfile(root, 'results', 'expert_multitarget_sweep');
end
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
addpath(scripts_dir);

source_imu_dir = 'C:\Users\LiuShuo\Desktop\results\IMU';
clutter_file = 'D:\Semester study\Owncode\ULA\single array element\s_clutter_full.mat';
manifest = prepare_real_imu(source_imu_dir, data_dir);
imu_file = manifest.replay_file{1};

target_speeds = [0.5, 1, 2, 3, 4, 6, 8, 10, 12, 15];
snr_list = [-25, -15, -5];
base_pos = {
    [180, 260, 90];
    [280, 180, 110];
    [520, 380, 85];
    [150, 200, 105];
    [400, 220, 95];
    [300, 180, 100];
    [300, 240, 100];
    [400, 100, 100];
    [400, 300, 100];
    [600, 400, 110]
};
n_target = numel(target_speeds);
targets = repmat(struct('coor0', [0, 0, 0], 'vt', [0, 0, 0], ...
    'A', [0, 0, 0], 'type', 'const'), 1, n_target);
for k = 1:n_target
    targets(k).coor0 = base_pos{k};
    dir_sign = (-1)^(k+1);
    targets(k).vt = [dir_sign * target_speeds(k), 0.2 * dir_sign * target_speeds(k), 0];
end

scheme_names = {'Ap', 'Bp', 'C', 'E', 'F', 'D', 'DnoIMU'};
scheme_labels = {
    'MTI';
    'Wide notch';
    'RD soft notch';
    'PCA subspace';
    'Hough ridge';
    'Proposed';
    'D-noIMU'
};
rows = table();

for si_sweep = 1:numel(snr_list)
    for ti_sweep = 1:n_target
        cfg = struct();
        cfg.batch_mode = true;
        cfg.save_figures = false;
        cfg.clutter_file = clutter_file;
        cfg.imu_file = imu_file;
        cfg.snr_db = snr_list(si_sweep);
        cfg.max_frames = 4;
        cfg.mid_frame = 2;
        cfg.random_seed = 3000 + 100 * si_sweep + ti_sweep;
        cfg.pfa_list = logspace(-4, -1, 12);
        cfg.multi_targets = targets;
        cfg.target_index = ti_sweep;
        cfg.radar_offset = [0, 0, 5];

        fprintf('[sweep] SNR=%d dB, target=%d/%d, speed=%.1f m/s\n', ...
            snr_list(si_sweep), ti_sweep, n_target, target_speeds(ti_sweep));
        EXPERIMENT_CFG = cfg; %#ok<NASGU>
        run(fullfile(scripts_dir, 'ob_moco_v4_ablation_experiment.m'));

        for mi_sweep = 1:numel(scheme_names)
            sn = scheme_names{mi_sweep};
            row = table();
            row.snr_db = snr_list(si_sweep);
            row.target_index = ti_sweep;
            row.target_speed = target_speeds(ti_sweep);
            row.method = string(scheme_names{mi_sweep});
            row.label = string(scheme_labels{mi_sweep});
            row.SCR_dB = mean(schemes.(sn).SCR, 'omitnan');
            valid = peak_found & isfinite(schemes.(sn).TERI);
            row.TERI = mean(schemes.(sn).TERI(valid), 'omitnan');
            row.FAR = far_avg(mi_sweep);
            row.Pd = mean(schemes.(sn).Pd, 'omitnan');
            roc_far = mean(roc_results.(sn).FAR_vs_Pfa, 1, 'omitnan');
            roc_pd = mean(roc_results.(sn).Pd_vs_Pfa, 1, 'omitnan');
            valid_roc = isfinite(roc_far) & isfinite(roc_pd) & roc_far > 0;
            if nnz(valid_roc) >= 2
                row.AUC = trapz(log10(roc_far(valid_roc)), roc_pd(valid_roc)) / ...
                    (log10(max(roc_far(valid_roc))) - log10(min(roc_far(valid_roc))));
            else
                row.AUC = NaN;
            end
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end

out_csv = fullfile(results_dir, 'expert_multitarget_sweep_results.csv');
out_mat = fullfile(results_dir, 'expert_multitarget_sweep_results.mat');
writetable(rows, out_csv);
save(out_mat, 'rows', 'target_speeds', 'snr_list', 'targets', 'imu_file', '-v7.3');

plot_sweep_heatmaps(rows, target_speeds, snr_list, scheme_names, scheme_labels, results_dir);
fprintf('[sweep] Results saved to %s\n', out_csv);
end

function plot_sweep_heatmaps(rows, target_speeds, snr_list, scheme_names, scheme_labels, results_dir)
methods = {'C', 'E', 'DnoIMU', 'D'};
method_labels = {'C: RD soft notch', 'E: PCA subspace', ...
    'D-noIMU: TF ridge', 'D: Proposed'};
metrics = {'SCR_dB', 'TERI', 'FAR', 'Pd', 'AUC'};
metric_labels = {'SCR (dB)', 'TERI', 'FAR', 'P_d', 'AUC'};

for m_idx = 1:numel(metrics)
    metric = metrics{m_idx};
    fig = figure('Color', 'w', 'Position', [80, 80, 1500, 900]);
    for mi = 1:numel(methods)
        method = methods{mi};
        Z = nan(numel(target_speeds), numel(snr_list));
        for si = 1:numel(snr_list)
            for ti = 1:numel(target_speeds)
                mask = rows.snr_db == snr_list(si) & rows.target_index == ti & ...
                    rows.method == string(method);
                if any(mask)
                    Z(ti, si) = rows.(metric)(find(mask, 1));
                end
            end
        end
        subplot(2, 2, mi);
        imagesc(Z); colorbar; axis xy;
        set(gca, 'XTick', 1:numel(snr_list), 'XTickLabel', snr_list, ...
            'YTick', 1:numel(target_speeds), 'YTickLabel', target_speeds, ...
            'FontSize', 10);
        xlabel('SNR (dB)'); ylabel('Target speed (m/s)');
        title(method_labels{mi}, 'FontSize', 11);
    end
    sgtitle(['Multi-target sweep: ' metric_labels{m_idx}], 'FontSize', 13);
    png_file = fullfile(results_dir, ['sweep_' metric '.png']);
    fig_file = fullfile(results_dir, ['sweep_' metric '.fig']);
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
    set(fig, 'Color', 'w', 'PaperPositionMode', 'auto');
    print(fig, png_file, '-dpng', '-r300');
    saveas(fig, fig_file);
    close(fig);
end
end
