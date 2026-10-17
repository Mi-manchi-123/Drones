function out_csv = run_expert_required_experiment()
%RUN_EXPERT_REQUIRED_EXPERIMENT Real-IMU replay with PCA and no-IMU baselines.

this_dir = fileparts(mfilename('fullpath'));
if exist(fullfile(this_dir, 'ob_moco_v4_ablation_experiment.m'), 'file')
    scripts_dir = this_dir;
    root = this_dir;
    data_dir = fullfile(root, 'imu_replay');
    results_dir = fullfile(root, 'results_expert_required_real_imu');
else
    scripts_dir = this_dir;
    root = fileparts(this_dir);
    data_dir = fullfile(root, 'data', 'imu_replay');
    results_dir = fullfile(root, 'results', 'expert_required_real_imu');
end
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
addpath(scripts_dir);

source_imu_dir = 'C:\Users\LiuShuo\Desktop\results\IMU';
clutter_file = 'D:\Semester study\Owncode\ULA\single array element\s_clutter_full.mat';
manifest = prepare_real_imu(source_imu_dir, data_dir);
imu_file = manifest.replay_file{1};

cfg = struct();
cfg.batch_mode = true;
cfg.save_figures = false;
cfg.clutter_file = clutter_file;
cfg.imu_file = imu_file;
cfg.snr_db = -15;
cfg.target_vt = [3.2, 2.4, 0];
cfg.radar_offset = [0, 0, 5];
cfg.random_seed = 2026;
cfg.pfa_list = logspace(-5, -1, 20);
cfg.max_frames = 24;
cfg.mid_frame = 12;

EXPERIMENT_CFG = cfg; %#ok<NASGU>
run(fullfile(scripts_dir, 'ob_moco_v4_ablation_experiment.m'));

scheme_names = {'Ap', 'Bp', 'C', 'E', 'F', 'D', 'DnoIMU'};
scheme_labels = {
    'A: MTI', ...
    'B: Wide notch', ...
    'C: RD soft notch', ...
    'E: PCA subspace', ...
    'F: Hough ridge', ...
    'D: Proposed IMU+Viterbi', ...
    'D-noIMU: TF ridge without IMU'
};

rows = table();
for i = 1:numel(scheme_names)
    sn = scheme_names{i};
    row = table();
    row.method = string(scheme_names{i});
    row.label = string(scheme_labels{i});
    row.SCR_dB = mean(schemes.(sn).SCR, 'omitnan');
    valid = peak_found & isfinite(schemes.(sn).TERI);
    row.TERI = mean(schemes.(sn).TERI(valid), 'omitnan');
    row.FAR = far_avg(i);
    row.Pd = mean(schemes.(sn).Pd, 'omitnan');
    roc_far = mean(roc_results.(sn).FAR_vs_Pfa, 1, 'omitnan');
    roc_pd  = mean(roc_results.(sn).Pd_vs_Pfa, 1, 'omitnan');
    valid_roc = isfinite(roc_far) & isfinite(roc_pd) & roc_far > 0;
    if nnz(valid_roc) >= 2
        row.AUC = trapz(log10(roc_far(valid_roc)), roc_pd(valid_roc)) / ...
            (log10(max(roc_far(valid_roc))) - log10(min(roc_far(valid_roc))));
    else
        row.AUC = NaN;
    end
    rows = [rows; row]; %#ok<AGROW>
end

out_csv = fullfile(results_dir, 'expert_required_real_imu_results.csv');
out_mat = fullfile(results_dir, 'expert_required_real_imu_results.mat');
writetable(rows, out_csv);
save(out_mat, 'rows', 'schemes', 'roc_results', 'PFA_LIST', 'cfg', 'MID', '-v7.3');

% ---- Processed RD maps for the two newly added methods ----
if exist('MID', 'var') && isfield(MID, 'rd_E') && isfield(MID, 'rd_Dno')
    rd_items = {
        MID.rd_E,  'PCA subspace processed RD',   'rd_pca_subspace.png';
        MID.rd_Dno,'No-IMU TF ridge processed RD','rd_tf_ridge_noimu.png'
    };
    for k = 1:size(rd_items, 1)
        f = figure('Color', 'w', 'Position', [100, 100, 900, 650]);
        imagesc(range_axis, v_tf, 20*log10(abs(rd_items{k,1})+eps));
        axis xy; colorbar; clim([-10, 100]); hold on;
        if isfield(MID, 'r_ref') && isfield(MID, 'd_ref')
            plot(range_axis(MID.r_ref), v_tf(MID.d_ref), 'r^', ...
                'MarkerSize', 8, 'LineWidth', 1.4);
        end
        xlabel('Range (m)'); ylabel('Speed (m/s)');
        title(rd_items{k, 2});
        set(findall(f, '-property', 'FontName'), 'FontName', 'Arial');
        set(gca, 'FontSize', 11, 'LineWidth', 0.9);
        grid on;
        print(f, fullfile(results_dir, rd_items{k, 3}), '-dpng', '-r300');
        saveas(f, fullfile(results_dir, strrep(rd_items{k, 3}, '.png', '.fig')));
        close(f);
    end
end

plot_methods = {'C', 'E', 'DnoIMU', 'D'};
plot_labels = {'RD soft notch', 'PCA subspace', 'No-IMU ridge', 'Proposed'};
idx = zeros(1, numel(plot_methods));
for i = 1:numel(plot_methods)
    idx(i) = find(strcmp(scheme_names, plot_methods{i}), 1);
end

fig = figure('Color', 'w', 'Position', [100, 100, 1200, 720]);
subplot(2,2,1);
bar(rows.SCR_dB(idx), 'FaceColor', [0.25 0.55 0.85]); grid on;
set(gca, 'XTick', 1:numel(idx), 'XTickLabel', plot_labels, 'FontSize', 11);
ylabel('SCR (dB)'); title('SCR');

subplot(2,2,2);
bar(rows.TERI(idx), 'FaceColor', [0.85 0.45 0.20]); grid on;
set(gca, 'XTick', 1:numel(idx), 'XTickLabel', plot_labels, 'FontSize', 11);
ylabel('TERI'); title('Target energy retention');

subplot(2,2,3);
bar(rows.Pd(idx), 'FaceColor', [0.30 0.65 0.35]); grid on; ylim([0, 1.05]);
set(gca, 'XTick', 1:numel(idx), 'XTickLabel', plot_labels, 'FontSize', 11);
ylabel('P_d'); title('Detection probability');

subplot(2,2,4);
bar(rows.FAR(idx), 'FaceColor', [0.75 0.25 0.25]); grid on;
set(gca, 'XTick', 1:numel(idx), 'XTickLabel', plot_labels, 'FontSize', 11);
ylabel('FAR'); title('False alarm rate');

sgtitle('Real-IMU replay: expert-required comparison', 'FontSize', 13);
print(fig, fullfile(results_dir, 'expert_required_comparison.png'), '-dpng', '-r300');
saveas(fig, fullfile(results_dir, 'expert_required_comparison.fig'));
close(fig);

fprintf('[expert] Results saved to %s\n', out_csv);
end
