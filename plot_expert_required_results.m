function plot_expert_required_results(results_dir)
%PLOT_EXPERT_REQUIRED_RESULTS Plot all metrics and RD comparisons.

if nargin < 1 || isempty(results_dir)
    root = fileparts(mfilename('fullpath'));
    root = fileparts(root);
    results_dir = fullfile(root, 'results', 'expert_required_real_imu');
end

mat_file = fullfile(results_dir, 'expert_required_real_imu_results.mat');
if ~exist(mat_file, 'file')
    error('Result MAT file not found: %s', mat_file);
end
D = load(mat_file, 'rows', 'schemes', 'roc_results', 'PFA_LIST', 'MID', 'cfg');
rows = D.rows;
schemes = D.schemes;
roc_results = D.roc_results;
PFA_LIST = D.PFA_LIST;
MID = D.MID;

labels = rows.label;
short_labels = {'A MTI', 'B Wide notch', 'C RD soft', 'E PCA', ...
    'F Hough', 'D Proposed', 'D-noIMU'};
metric_specs = {
    'SCR_dB', 'SCR (dB)', 'SCR';
    'TERI', 'TERI', 'Target energy retention';
    'FAR', 'FAR', 'False alarm rate';
    'Pd', 'P_d', 'Detection probability';
    'AUC', 'AUC', 'ROC area'
};

% Individual metric figures.
for i = 1:size(metric_specs, 1)
    key = metric_specs{i, 1};
    fig = figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    bar(rows.(key), 'FaceColor', [0.25, 0.55, 0.85]); grid on;
    set(gca, 'XTick', 1:numel(short_labels), 'XTickLabel', short_labels, ...
        'XTickLabelRotation', 25, 'FontSize', 11);
    ylabel(metric_specs{i, 2});
    title(metric_specs{i, 3}, 'FontSize', 12);
    export_figure(fig, fullfile(results_dir, ['expert_required_metric_' key '.png']), ...
        fullfile(results_dir, ['expert_required_metric_' key '.fig']));
    close(fig);
end

% Combined metric panel.
fig = figure('Color', 'w', 'Position', [80, 80, 1500, 900]);
for i = 1:size(metric_specs, 1)
    subplot(2, 3, i);
    bar(rows.(metric_specs{i, 1}), 'FaceColor', [0.25, 0.55, 0.85]); grid on;
    set(gca, 'XTick', 1:numel(short_labels), 'XTickLabel', short_labels, ...
        'XTickLabelRotation', 30, 'FontSize', 10);
    ylabel(metric_specs{i, 2});
    title(metric_specs{i, 3}, 'FontSize', 11);
end
subplot(2, 3, 6);
plot_roc_axes(gca, roc_results, PFA_LIST);
sgtitle('Expert-required real-IMU experiment: all evaluation metrics', ...
    'FontSize', 13);
export_figure(fig, fullfile(results_dir, 'expert_required_all_metrics.png'), ...
    fullfile(results_dir, 'expert_required_all_metrics.fig'));
close(fig);

% ROC curves.
fig = figure('Color', 'w', 'Position', [100, 100, 900, 650]);
plot_roc_axes(gca, roc_results, PFA_LIST);
title('ROC curves from recorded-IMU replay', 'FontSize', 12);
export_figure(fig, fullfile(results_dir, 'expert_required_roc.png'), ...
    fullfile(results_dir, 'expert_required_roc.fig'));
close(fig);

% Processed RD comparison.
if ~isempty(fieldnames(MID))
    fig = figure('Color', 'w', 'Position', [50, 50, 1500, 900]);
    rd_data = {MID.rd_pre, MID.rd_C, MID.rd_E, MID.rd_Dno, MID.rd_D};
    rd_titles = {'Input RD', 'C: RD soft notch', 'E: PCA subspace', ...
        'D-noIMU: TF ridge', 'D: Proposed'};
    for i = 1:numel(rd_data)
        subplot(2, 3, i);
        imagesc(20*log10(abs(rd_data{i})+eps)); axis xy; colorbar;
        clim([-10, 100]); grid on;
        set(gca, 'FontSize', 10);
        xlabel('Range bin'); ylabel('Doppler bin');
        title(rd_titles{i}, 'FontSize', 11);
    end
    sgtitle('Processed RD comparison: expert-required baselines', 'FontSize', 13);
    export_figure(fig, fullfile(results_dir, 'expert_required_rd_comparison.png'), ...
        fullfile(results_dir, 'expert_required_rd_comparison.fig'));
    close(fig);
end

fprintf('[plot_expert_required_results] Figures saved in %s\n', results_dir);
end

function plot_roc_axes(ax, roc_results, pfa_list)
scheme_names = {'Ap', 'Bp', 'C', 'E', 'F', 'D', 'DnoIMU'};
labels = {'A MTI', 'B Wide notch', 'C RD soft', 'E PCA', ...
    'F Hough', 'D Proposed', 'D-noIMU'};
colors = lines(numel(scheme_names));
hold(ax, 'on');
for i = 1:numel(scheme_names)
    sn = scheme_names{i};
    far = mean(roc_results.(sn).FAR_vs_Pfa, 1, 'omitnan');
    pd = mean(roc_results.(sn).Pd_vs_Pfa, 1, 'omitnan');
    valid = isfinite(far) & isfinite(pd) & far > 0;
    if any(valid)
        plot(ax, far(valid), pd(valid), '-o', 'Color', colors(i, :), ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', labels{i});
    end
end
set(ax, 'XScale', 'log', 'FontSize', 10);
xlabel(ax, 'Actual FAR'); ylabel(ax, 'Detection probability P_d');
legend(ax, 'Location', 'southeast', 'FontSize', 8);
grid(ax, 'on'); xlim(ax, [1e-6, 1]); ylim(ax, [0, 1.05]);
end

function export_figure(fig, png_file, fig_file)
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(fig, 'Color', 'w', 'PaperPositionMode', 'auto');
print(fig, png_file, '-dpng', '-r300');
saveas(fig, fig_file);
end
