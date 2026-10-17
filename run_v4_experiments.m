%% =========================================================================
%  run_v4_experiments.m
%  多场景消融实验运行脚本 — 配合 ob_moco_v4_ablation.m
%
%  功能：
%    1. 多 SNR 场景（-20, -15, -10 dB）
%    2. 多目标速度（低速 2 m/s、中速 4 m/s、高速 6 m/s）
%    3. 蒙特卡洛多次运行（可选）
%    4. 输出结构化结果到 MAT 文件供后续分析
%
%  用法：直接运行（设置对应参数后自动调用 v4 核心处理函数）
% =========================================================================
clear; clc; close all;
fprintf('============================================================\n');
fprintf('  多场景消融实验运行脚本\n');
fprintf('============================================================\n\n');

%% =========================================================================
% 实验配置
% =========================================================================

% SNR 列表 (dB)
SNR_LIST = [-20, -15, -10];

% 目标速度列表（径向速度近似值，m/s）
% 每个元素为 [vx, vy, vz] 格式
TGT_VT_LIST = [
    1.5, 1.3, 0;   % 低速 ~2 m/s
    3.2, 2.4, 0;   % 中速 ~4 m/s（原 v3 默认）
    5.0, 3.3, 0    % 高速 ~6 m/s
];

% 是否做蒙特卡洛（多次独立随机杂波和噪声实现）
N_MONTE_CARLO = 1;    % 设为 5-10 做统计均值，设为 1 做单次

% 方案名称（与 v4 一致）
SCHEME_NAMES = {'Ap', 'Bp', 'C', 'E', 'F', 'D', 'DnoIMU'};
N_SCHEMES = length(SCHEME_NAMES);
N_SNR = length(SNR_LIST);
N_SPD = size(TGT_VT_LIST, 1);

%% =========================================================================
% 结果存储
% =========================================================================

% 主结构：results.(scheme_lower){snr_idx, spd_idx, mc_idx}.SCR ...
% 汇总结构：results.summary.SCR_mean(snr_idx, spd_idx, scheme_idx) ...
results = struct();
results.cfg.SNR_LIST  = SNR_LIST;
results.cfg.TGT_VT_LIST = TGT_VT_LIST;
results.cfg.N_MONTE_CARLO = N_MONTE_CARLO;

% 初始化方案结果单元数组
for sc = 1:N_SCHEMES
    sn = lower(SCHEME_NAMES{sc});
    results.(sn) = cell(N_SNR, N_SPD, N_MONTE_CARLO);
end
results.peak_found = cell(N_SNR, N_SPD, N_MONTE_CARLO);
results.mu_c = cell(N_SNR, N_SPD, N_MONTE_CARLO);
results.sigma_c = cell(N_SNR, N_SPD, N_MONTE_CARLO);
results.errors = {};

fprintf('实验配置:\n');
fprintf('  SNR 列表: ');
fprintf('%.0f ', SNR_LIST); fprintf('dB\n');
fprintf('  目标速度: ');
for s = 1:N_SPD
    fprintf('(%.1f,%.1f,%.1f) ', TGT_VT_LIST(s,1), TGT_VT_LIST(s,2), TGT_VT_LIST(s,3));
end
fprintf('\n');
fprintf('  蒙特卡洛次数: %d\n', N_MONTE_CARLO);
fprintf('  方案数: %d\n\n', N_SCHEMES);

%% =========================================================================
% 实验主循环
% =========================================================================
total_runs = N_SNR * N_SPD * N_MONTE_CARLO;
run_count = 0;

for si = 1:N_SNR
    SNR_db = SNR_LIST(si);
    fprintf('━━━ SNR = %.0f dB ━━━\n', SNR_db);

    for vi = 1:N_SPD
        tgt_vt = TGT_VT_LIST(vi, :);
        fprintf('  ━ 目标速度 (%.1f,%.1f,%.1f) m/s ━\n', tgt_vt(1), tgt_vt(2), tgt_vt(3));

        for mi = 1:N_MONTE_CARLO
            run_count = run_count + 1;
            fprintf('    蒙特卡洛 #%d/%d ... ', mi, N_MONTE_CARLO);

            % ── 调用 v4 核心 —— 注意：
            %   ob_moco_v4_ablation.m 是独立脚本，不能作为函数调用。
            %   此处使用两种策略：
            %     (1) 直接修改全局变量后 run(script) — 需要 v4 支持参数化
            %     (2) 重构 v4 为函数形式 — 更干净但改动大
            %
            %  当前采用策略 (1)：通过 assignin/evalin 注入参数后 run
            % ──

            % 设置目标速度
            assignin('base', 'TGT_VT_OPTIONS', tgt_vt);
            assignin('base', 'SNR_db_input', SNR_db);

            % 运行核心脚本（捕获输出到 diary）
            diary('v4_run_log.txt');
            try
                ob_moco_v4_ablation;
                % 收集结果（v4 脚本在 base workspace 中留下变量）
                SCR_data  = evalin('base', 'schemes');
                ROC_data  = evalin('base', 'roc_results');
                peak_ok   = evalin('base', 'peak_found');
                mu_data   = evalin('base', 'mu_c_all');
                sigma_data= evalin('base', 'sigma_c_all');

                % 保存到结果结构
                for sc = 1:N_SCHEMES
                    sn = SCHEME_NAMES{sc};
                    sn_lower = lower(sn);
                    results.(sn_lower){si, vi, mi}.SCR   = SCR_data.(sn).SCR;
                    results.(sn_lower){si, vi, mi}.TERI  = SCR_data.(sn).TERI;
                    results.(sn_lower){si, vi, mi}.BW3   = SCR_data.(sn).BW3;
                    results.(sn_lower){si, vi, mi}.pts   = SCR_data.(sn).pts;
                    results.(sn_lower){si, vi, mi}.Pd    = SCR_data.(sn).Pd;
                    results.(sn_lower){si, vi, mi}.ROC_Pd   = ROC_data.(sn).Pd_vs_Pfa;
                    results.(sn_lower){si, vi, mi}.ROC_FAR  = ROC_data.(sn).FAR_vs_Pfa;
                end
                results.peak_found{si, vi, mi} = peak_ok;
                results.mu_c{si, vi, mi} = mu_data;
                results.sigma_c{si, vi, mi} = sigma_data;

                fprintf('完成 ✓\n');
            catch ME
                fprintf('失败 ✗: %s\n', ME.message);
                % 保存部分结果（防御：被调脚本可能 clear 了工作空间）
                if ~exist('results', 'var')
                    results = struct();
                    results.cfg.SNR_LIST = SNR_LIST;
                    results.cfg.TGT_VT_LIST = TGT_VT_LIST;
                    results.cfg.N_MONTE_CARLO = N_MONTE_CARLO;
                end
                if ~isfield(results, 'errors')
                    results.errors = {};
                end
                results.errors{end+1} = struct('snr_idx',si,'spd_idx',vi,'mc_idx',mi,'msg',ME.message);
            end
            diary('off');
        end
    end
end

fprintf('\n实验完成。总运行次数: %d\n', run_count);

%% =========================================================================
% 结果汇总统计
% =========================================================================
fprintf('\n============================================================\n');
fprintf('  结果汇总\n');
fprintf('============================================================\n');

% 对每个 (SNR, 速度, 方案) 计算均值
summary = struct();

for si = 1:N_SNR
    for vi = 1:N_SPD
        % 合并所有蒙特卡洛运行
        for sc = 1:N_SCHEMES
            sn = SCHEME_NAMES{sc};
            sn_lower = lower(sn);

            % 收集所有 MC 的 SCR 均值
            scr_mc = zeros(N_MONTE_CARLO, 1);
            teri_mc = zeros(N_MONTE_CARLO, 1);
            far_mc = zeros(N_MONTE_CARLO, 1);

            for mi = 1:N_MONTE_CARLO
                if ~isempty(results.(sn_lower){si, vi, mi})
                    scr_mc(mi)  = mean(results.(sn_lower){si, vi, mi}.SCR, 'omitnan');
                    teri_mc(mi) = mean(results.(sn_lower){si, vi, mi}.TERI, 'omitnan');
                    % FAR: 总检测点数 / (N_slow * N_fast)
                    far_mc(mi)  = mean(results.(sn_lower){si, vi, mi}.pts) / (256*1024);
                end
            end

            summary.SCR_mean(si, vi, sc)  = mean(scr_mc, 'omitnan');
            summary.SCR_std(si, vi, sc)   = std(scr_mc, 'omitnan');
            summary.TERI_mean(si, vi, sc) = mean(teri_mc, 'omitnan');
            summary.FAR_mean(si, vi, sc)  = mean(far_mc, 'omitnan');
        end
    end
end

% ── 打印汇总表 ──
fprintf('\n──────────────────────────────────────────────────────\n');
fprintf('  SCR 均值 (dB) — 行=SNR, 列=方案, 子表=速度\n');
fprintf('──────────────────────────────────────────────────────\n');
for vi = 1:N_SPD
    fprintf('\n  目标速度 #%d: (%.1f,%.1f,%.1f) m/s\n', vi, ...
        TGT_VT_LIST(vi,1), TGT_VT_LIST(vi,2), TGT_VT_LIST(vi,3));
    fprintf('  %-8s', 'SNR(dB)');
    for sc = 1:N_SCHEMES; fprintf(' %8s', SCHEME_NAMES{sc}); end
    fprintf('\n');
    for si = 1:N_SNR
        fprintf('  %-8.0f', SNR_LIST(si));
        for sc = 1:N_SCHEMES
            fprintf(' %8.2f', summary.SCR_mean(si, vi, sc));
        end
        fprintf('\n');
    end
end

% ── 关键消融结论表 ──
fprintf('\n──────────────────────────────────────────────────────\n');
fprintf('  关键交叉对比（SCR 差异 ΔdB）\n');
fprintf('──────────────────────────────────────────────────────\n');
for si = 1:N_SNR
    fprintf('\n  SNR = %.0f dB:\n', SNR_LIST(si));
    fprintf('  %-28s %8s %8s %8s\n', '对比项', '速度#1', '速度#2', '速度#3');
    % (1) 时频面 vs RD面: C vs D-noIMU
    diff_tf_vs_rd = squeeze(summary.SCR_mean(si, :, 7) - summary.SCR_mean(si, :, 3));
    fprintf('  %-28s %8.2f %8.2f %8.2f\n', '(1) 时频面vsRD面(C↔Dno)', diff_tf_vs_rd);
    % (2) Viterbi vs Hough: F vs D-noIMU
    diff_vit_vs_hough = squeeze(summary.SCR_mean(si, :, 7) - summary.SCR_mean(si, :, 5));
    fprintf('  %-28s %8.2f %8.2f %8.2f\n', '(2) Viterbi vs Hough(F↔Dno)', diff_vit_vs_hough);
    % (3) IMU先验: D-noIMU vs D
    diff_imu = squeeze(summary.SCR_mean(si, :, 6) - summary.SCR_mean(si, :, 7));
    fprintf('  %-28s %8.2f %8.2f %8.2f\n', '(3) IMU先验贡献(Dno↔D)', diff_imu);
end

%% =========================================================================
% 保存结果
% =========================================================================
save('v4_experiment_results.mat', 'results', 'summary', '-v7.3');
fprintf('\n结果已保存到 v4_experiment_results.mat\n');

%% =========================================================================
% 图：多场景消融汇总
% =========================================================================
figure('Name','多场景消融汇总','Position',[30,30,1400,800]);

% 子图1：SNR vs SCR（中速场景）
subplot(2,3,1);
vi_mid = 2;  % 中速
for sc = 1:N_SCHEMES
    plot(SNR_LIST, summary.SCR_mean(:, vi_mid, sc), 'o-', ...
        'LineWidth', 1.2+0.3*(sc>=6), 'DisplayName', SCHEME_NAMES{sc}); hold on;
end
xlabel('SNR (dB)'); ylabel('SCR (dB)');
title(sprintf('SCR vs SNR (速度 #%d)', vi_mid));
legend('Location','best','FontSize',7); grid on;

% 子图2：目标速度 vs SCR（SNR=-15dB）
subplot(2,3,2);
si_mid = 2;  % SNR = -15 dB
speed_labels = {'低速','中速','高速'};
for sc = 1:N_SCHEMES
    plot(1:N_SPD, squeeze(summary.SCR_mean(si_mid, :, sc)), 'o-', ...
        'LineWidth', 1.2+0.3*(sc>=6), 'DisplayName', SCHEME_NAMES{sc}); hold on;
end
set(gca, 'XTick', 1:N_SPD, 'XTickLabel', speed_labels);
xlabel('目标速度'); ylabel('SCR (dB)');
title(sprintf('SCR vs 目标速度 (SNR=%.0f dB)', SNR_LIST(si_mid)));
legend('Location','best','FontSize',7); grid on;

% 子图3：消融对比柱状图（SNR=-15dB, 中速）
subplot(2,3,3);
bar_data = squeeze(summary.SCR_mean(si_mid, vi_mid, :));
bar(1:N_SCHEMES, bar_data, 0.6, 'FaceColor', [0.3 0.6 0.9]); hold on;
for sc = 1:N_SCHEMES
    text(sc, bar_data(sc)+0.3, sprintf('%.1f', bar_data(sc)), ...
        'Hor','center','Vert','bottom','FontSize',8);
end
set(gca,'XTick',1:N_SCHEMES,'XTickLabel',SCHEME_NAMES,'XTickLabelRotation',30);
ylabel('SCR (dB)'); title(sprintf('消融对比 (SNR=%.0f, 中速)', SNR_LIST(si_mid)));
grid on;
annotation('textarrow',[0.22,0.35],[0.5,0.5],'String','RD域链','FontSize',8);
annotation('textarrow',[0.52,0.75],[0.5,0.5],'String','时频域链','FontSize',8);

% 子图4：TERI 对比
subplot(2,3,4);
teri_data = squeeze(summary.TERI_mean(si_mid, vi_mid, 3:N_SCHEMES));
bar(1:length(teri_data), teri_data, 0.6, 'FaceColor', [0.9 0.5 0.2]); hold on;
yline(1.0, 'k--', 'LineWidth', 1.0);
for sc = 1:length(teri_data)
    text(sc, teri_data(sc)+0.02, sprintf('%.2f', teri_data(sc)), ...
        'Hor','center','Vert','bottom','FontSize',8);
end
set(gca,'XTick',1:length(teri_data),'XTickLabel',SCHEME_NAMES(3:N_SCHEMES),'XTickLabelRotation',30);
ylabel('TERI'); title('目标能量保持率'); grid on;

% 子图5：FAR 对比
subplot(2,3,5);
far_data = squeeze(summary.FAR_mean(si_mid, vi_mid, :));
bar(1:N_SCHEMES, far_data*1e4, 0.6, 'FaceColor', [0.7 0.2 0.2]); hold on;
set(gca,'XTick',1:N_SCHEMES,'XTickLabel',SCHEME_NAMES,'XTickLabelRotation',30);
ylabel('FAR (×10^{-4})'); title('虚警率对比'); grid on;

% 子图6：关键Δ表
subplot(2,3,6); axis off;
d1 = squeeze(summary.SCR_mean(si_mid, vi_mid, 7) - summary.SCR_mean(si_mid, vi_mid, 3));
d2 = squeeze(summary.SCR_mean(si_mid, vi_mid, 7) - summary.SCR_mean(si_mid, vi_mid, 5));
d3 = squeeze(summary.SCR_mean(si_mid, vi_mid, 6) - summary.SCR_mean(si_mid, vi_mid, 7));
text(0.1,0.9,sprintf([ ...
    '关键结论 (SNR=%.0f dB, 中速):\n\n', ...
    '(1) 时频面 vs RD面:\n', ...
    '    D-noIMU(%.1f) vs C(%.1f) = %.2f dB\n\n', ...
    '(2) Viterbi vs Hough:\n', ...
    '    D-noIMU(%.1f) vs F(%.1f) = %.2f dB\n\n', ...
    '(3) IMU先验贡献:\n', ...
    '    D(%.1f) vs D-noIMU(%.1f) = %+.2f dB\n'], ...
    SNR_LIST(si_mid), ...
    summary.SCR_mean(si_mid, vi_mid, 7), summary.SCR_mean(si_mid, vi_mid, 3), d1, ...
    summary.SCR_mean(si_mid, vi_mid, 7), summary.SCR_mean(si_mid, vi_mid, 5), d2, ...
    summary.SCR_mean(si_mid, vi_mid, 6), summary.SCR_mean(si_mid, vi_mid, 7), d3), ...
    'FontSize',11,'VerticalAlignment','top');

sgtitle('多场景消融实验汇总 — SNR×目标速度×蒙特卡洛', ...
    'FontSize',13,'FontWeight','bold');

saveas(gcf, 'v4_ablation_summary.png');
fprintf('汇总图已保存到 v4_ablation_summary.png\n');

fprintf('\n============================================================\n');
fprintf('  实验结束。运行 run_v4_analysis.m 进行深入分析。\n');
fprintf('============================================================\n');
