%% =========================================================================
%  ob_moco_v3_tf_ridge.m
%  单通道车载雷达：慢时间时频(STFT)动态杂波建模 + Viterbi 脊线追踪 +
%                   目标保持型时频掩模滤波
%
%  在 ob_moco_v2_dynamic_clutter.m 基础上完善：平台/目标参数与杂波、回波
%  信号生成完全沿用原代码（PART 0–6），处理链替换为论文所述的时频流水线。
%
%  四方案对比（统一到 N_slow 点 FFT 的同一条多普勒速度轴 v_tf 上）：
%    A：仅平动补偿 + 传统 MTI 滤波      （二脉冲对消，基线）
%    B：MOCO + 宽硬陷波                 （过抑制基线，损伤慢目标）
%    C：MOCO + RD 域自适应软陷波        （原 v2 方法：每 CPI 单一 μ_c/σ_c）
%    D：MOCO + 慢时间 STFT + Viterbi 脊线追踪 + 目标保持型时频掩模（本文/论文方法）
%
%  核心新增（相对 v2）：
%    (1) 慢时间 STFT，把 CPI 内的杂波脊线"蜿蜒"分辨到时频面 (创新点1·时频表征)
%    (2) IMU 软先验 + Viterbi 动态规划在时频面追踪脊线           (创新点2)
%    (3) 平顶+高斯滚降的时频掩模 + 候选目标全程保护 + WOLA-ISTFT 重建 (创新点3)
%
%  评价指标：SCR / SCR增益 / 目标能量保持率TERI / 多普勒谱宽压缩 / FAR
%
%  说明：A/B/C/D 共用同一 v_tf 轴，便于隔离 "RD 软陷波(C)" 与 "时频脊线(D)" 的差异。
% =========================================================================
clear; clc; close all;

fprintf('============================================================\n');
fprintf('  时频脊线追踪 + 目标保持型时频掩模 — 仿真验证平台 (v3)\n');
fprintf('============================================================\n\n');

%% =========================================================================
% PART 0：公共雷达参数（沿用 ob_moco_v2）
% =========================================================================
Ob_time  = 2.0;             % 总观测时间 (s)
fc       = 13.5e9;          % 中心频率 (Hz)
B        = 20e6;            % 信号带宽
Tc       = 5e-6;            % Chirp 脉宽
PRI      = 200e-6;          % 脉冲重复间隔
fs       = 60e6;            % 快时间采样率
c        = 3e8;             % 光速
lambda   = c / fc;          % 波长 ≈ 0.0222 m
K        = B / Tc;          % 调频斜率

N_fast   = 1024;            % 距离门数
N_slow   = 256;             % CPI 脉冲数
overlap  = 200;             % 帧间重叠（增大以提高帧数平滑Pd曲线，不增加观测时长）
Num_chirp = floor(Ob_time / PRI);   % 总脉冲数 = 10000

% 距离轴
t_fast      = (0:N_fast-1) / fs;
range_axis  = t_fast * c / 2;
range_res   = c / (fs * 2);

% 速度轴（CPI 级 / FFT 后 fftshift 约定）—— 四方案统一使用 v_tf
step       = N_slow - overlap;
num_frames = floor((Num_chirp - N_slow) / step) + 1;
v_tf       = ((0:N_slow-1) - floor(N_slow/2)) * (1/PRI) / N_slow * lambda/2;   % m/s
v_tf       = v_tf(:).';

% 窗函数
win_range   = hanning(N_fast).';
win_doppler = hanning(N_slow);

fprintf('[参数] λ=%.4f m | 距离分辨率=%.2f m | 多普勒速度分辨率=%.3f m/s\n', ...
    lambda, range_res, v_tf(2)-v_tf(1));
fprintf('[参数] 总脉冲数=%d | CPI帧数=%d | 每帧 N_slow=%d | v_tf∈[%.1f, %.1f] m/s\n\n', ...
    Num_chirp, num_frames, N_slow, v_tf(1), v_tf(end));

%% =========================================================================
% PART 1：平台运动状态生成（ESKF数据 或 仿真替代）
% =========================================================================
USE_SIM_IMU = ~exist('Data_deduced_in_reverse.csv', 'file');

if ~USE_SIM_IMU
    fprintf('[数据] 加载 ESKF 反推文件...\n');
    eskf_data = readtable('Data_deduced_in_reverse.csv');
else
    fprintf('[仿真] 未找到ESKF文件，使用仿真平台运动模型生成姿态序列\n');
    t_sim    = (0:Num_chirp-1).' * PRI;
    v_plat   = 10;          % 前进速度 (m/s)
    heading  = 30;          % 初始航向角 (°)
    P_x = v_plat * cos(deg2rad(heading)) * t_sim;
    P_y = v_plat * sin(deg2rad(heading)) * t_sim;
    P_z = zeros(Num_chirp, 1);
    A_x = 3.0 * sin(2*pi*1.5*t_sim) + 0.2*sin(2*pi*25*t_sim + 0.7);  % 横滚 (°)
    A_y = 1.5 * sin(2*pi*1.2*t_sim + 1.0) + 0.1*sin(2*pi*30*t_sim);  % 俯仰 (°)
    A_z = heading + 0.5*sin(2*pi*0.5*t_sim);                          % 偏航 (°)
    eskf_data.t   = t_sim;
    eskf_data.P_x = P_x; eskf_data.P_y = P_y; eskf_data.P_z = P_z;
    eskf_data.A_x = A_x; eskf_data.A_y = A_y; eskf_data.A_z = A_z;
    fprintf('[仿真] 平台速度=%.1f m/s，横滚=3°(1.5Hz)+0.2°(25Hz)\n\n', v_plat);
end

%% =========================================================================
% PART 2：多目标参数（蜂群无人机模拟）
% =========================================================================
% 设计要求（为突出方案D对慢速目标的保持能力，增加慢速目标数量）：
%   (1) 5个慢速目标(2.5~3 m/s)：空间打散且位于杂波带边缘，用于验证方案D的保持优势
%   (2) 2个中速目标(4~5 m/s)：方向相反的对比对
%   (3) 3个快速目标(7~9 m/s)：远离杂波带，验证快速目标检测不受影响
%   (4) 慢速目标分布 x∈[150,550], y∈[180,400], 充分分散
% =========================================================================
N_tgt = 10;

% ── 慢速目标（2.5~3 m/s，空间高度分散）─────────────────────────
%  这些目标的径向速度接近杂波带边缘(|v_rad|≈1~3 m/s)，
%  用于验证方案D vs C 对慢速目标的保持差异

% 目标1：左前慢速 speed≈2.83 m/s
tgt(1).coor0 = [180, 260, 90];
tgt(1).vt    = [2.5, 1.2, 0.4];
tgt(1).A     = [0, 0, 0];
tgt(1).type  = 'const';

% 目标2：右下慢速 speed≈2.83 m/s（与T3方向大致相反）
tgt(2).coor0 = [280, 180, 110];
tgt(2).vt    = [-2.4, 1.5, 0.0];
tgt(2).A     = [0, 0, 0];
tgt(2).type  = 'const';

% 目标3：右上远距慢速 speed≈2.77 m/s
tgt(3).coor0 = [520, 380, 85];
tgt(3).vt    = [2.6, -0.9, 0.3];
tgt(3).A     = [0, 0, 0];
tgt(3).type  = 'const';

% 目标4：左下慢速 speed≈2.75 m/s
tgt(4).coor0 = [150, 200, 105];
tgt(4).vt    = [0.9, 2.6, 0.0];
tgt(4).A     = [0, 0, 0];
tgt(4).type  = 'const';

% 目标5：中前慢速(最接近杂波带) speed≈2.83 m/s
tgt(5).coor0 = [400, 220, 95];
tgt(5).vt    = [2.8, 0.4, 0.1];
tgt(5).A     = [0, 0, 0];
tgt(5).type  = 'const';

% ── 中速目标（4~5 m/s，对比用）──────────────────────────────────

% 目标6：正向中速 speed≈4.47 m/s
tgt(6).coor0 = [300, 180, 100];
tgt(6).vt    = [4.0, 2.0, 0.0];
tgt(6).A     = [0, 0, 0];
tgt(6).type  = 'const';

% 目标7：与T6反向 speed≈4.47 m/s
tgt(7).coor0 = [300, 240, 100];
tgt(7).vt    = [-4.0, -2.0, 0.0];
tgt(7).A     = [0, 0, 0];
tgt(7).type  = 'const';

% ── 快速目标（7~8 m/s，远离杂波带，验证检测不受影响）────────────

% 目标8：正x远离 speed=7.0 m/s
tgt(8).coor0 = [400, 100, 100];
tgt(8).vt    = [7.0, 0.0, 0.0];
tgt(8).A     = [0, 0, 0];
tgt(8).type  = 'const';

% 目标9：负x靠近 speed=7.0 m/s（与T8相反）
tgt(9).coor0 = [400, 300, 100];
tgt(9).vt    = [-7.0, 0.0, 0.0];
tgt(9).A     = [0, 0, 0];
tgt(9).type  = 'const';

% 目标10：斜向快速 speed≈7.81 m/s
tgt(10).coor0 = [600, 400, 110];
tgt(10).vt    = [6.0, 5.0, 0.0];
tgt(10).A     = [0, 0, 0];
tgt(10).type  = 'const';

% 速度约束验证
speeds = arrayfun(@(k) norm(tgt(k).vt), 1:N_tgt);
fprintf('[目标] 共 %d 个蜂群目标 | 速度范围 [%.2f, %.2f] m/s\n', ...
    N_tgt, min(speeds), max(speeds));
fprintf('        空间分布 x∈[%.0f,%.0f] y∈[%.0f,%.0f] z∈[%.0f,%.0f]\n', ...
    min(reshape([tgt.coor0],3,[]),[],2), max(reshape([tgt.coor0],3,[]),[],2));
fprintf('        慢速(T1~T5): %.2f/%.2f/%.2f/%.2f/%.2f m/s (杂波带边缘)\n', ...
    speeds(1), speeds(2), speeds(3), speeds(4), speeds(5));
fprintf('        中速(T6~T7): %.2f/%.2f m/s (方向相反对比对)\n', ...
    speeds(6), speeds(7));
fprintf('        快速(T8~T10): %.2f/%.2f/%.2f m/s (远离杂波带)\n\n', ...
    speeds(8), speeds(9), speeds(10));

radar_offset = [0, 0, 5];       % 雷达阵面相对车体质心的安装偏移 (m)

% 蜂群几何中心（用作 MOCO 参考点 / 杂波脊线建模参考）
ref_coor0 = mean(reshape([tgt.coor0], 3, []), 2)';
IN_coor0 = ref_coor0;   % 保持下游代码变量名兼容

%% =========================================================================
% PART 3：批量预插值，计算天线真实位置与速度
% =========================================================================
t_all   = (0:Num_chirp-1).' * PRI;
pos_all = interp1(eskf_data.t, [eskf_data.P_x, eskf_data.P_y, eskf_data.P_z], t_all, 'spline', 'extrap');
ang_all = interp1(eskf_data.t, [eskf_data.A_x, eskf_data.A_y, eskf_data.A_z], t_all, 'spline', 'extrap');

dt_fd   = 1e-3;
pos_nxt = interp1(eskf_data.t, [eskf_data.P_x, eskf_data.P_y, eskf_data.P_z], t_all + dt_fd, 'spline', 'extrap');
ang_nxt = interp1(eskf_data.t, [eskf_data.A_x, eskf_data.A_y, eskf_data.A_z], t_all + dt_fd, 'spline', 'extrap');

vel_all      = (pos_nxt - pos_all) / dt_fd;
ang_rate_all = deg2rad(ang_nxt - ang_all) / dt_fd;

ant_vel_all = zeros(Num_chirp, 3);
ant_pos_all = zeros(Num_chirp, 3);
R_all       = zeros(3, 3, Num_chirp);

for i = 1:Num_chirp
    phi   = deg2rad(ang_all(i,1));
    theta = deg2rad(ang_all(i,2));
    psi   = deg2rad(ang_all(i,3));
    R_i   = build_R_zyx_local(phi, theta, psi);
    R_all(:,:,i) = R_i;

    arm_inertial     = R_i * radar_offset.';
    ant_pos_all(i,:) = pos_all(i,:) + arm_inertial.';

    T_e2b = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
        0, cos(phi),            -sin(phi);
        0, sin(phi)/cos(theta),  cos(phi)/cos(theta)];

    omega_body     = T_e2b \ ang_rate_all(i,:).';
    omega_inertial = R_i * omega_body;
    v_arm          = cross(omega_inertial, arm_inertial);
    ant_vel_all(i,:) = vel_all(i,:) + v_arm.';
end

ant_ideal_all = pos_all + radar_offset;

%% =========================================================================
% PART 4：多目标回波生成（蜂群）
% =========================================================================
tgt_R0   = norm(ref_coor0 - ant_pos_all(1,:));   % 参考距离（蜂群中心）

s_rx_raw    = zeros(Num_chirp, N_fast);
target_pos  = zeros(Num_chirp, 3, N_tgt);        % 所有时×所有目标位置
dist_target_arr = zeros(Num_chirp, N_tgt);        % 所有时×所有目标距离

for i = 1:Num_chirp
    t_curr       = (i-1) * PRI;
    ant_pos      = ant_pos_all(i,:);
    for k = 1:N_tgt
        tgt_pos_k    = tgt(k).coor0 + t_curr * tgt(k).vt + 0.5 * tgt(k).A * t_curr^2;
        target_pos(i,:,k) = tgt_pos_k;
        dist_t       = norm(tgt_pos_k - ant_pos);
        dist_target_arr(i,k) = dist_t;

        range_amp = (tgt_R0 / dist_t)^2;
        tau       = 2 * dist_t / c;
        s_rx_raw(i,:) = s_rx_raw(i,:) + ...
            range_amp * rectpuls(t_fast - Tc/2 - tau, Tc) .* ...
            exp(1j*pi*K*(t_fast - tau).^2) .* exp(-1j*2*pi*fc*tau);
    end
end

%% =========================================================================
% PART 5：生成或加载杂波信号
% =========================================================================
USE_SIM_CLUTTER = ~exist('s_clutter_full.mat', 'file');

if ~USE_SIM_CLUTTER
    fprintf('[数据] 加载预生成杂波数据...\n');
    data           = load('s_clutter_full.mat', 's_clutter_full');
    s_clutter_full = data.s_clutter_full;
else
    fprintf('[仿真] 生成动态地面杂波（散射单元模型）...\n');
    s_clutter_full = sim_dynamic_clutter(Num_chirp, N_fast, fs, fc, PRI, K, Tc, ...
        t_fast, ant_pos_all, ant_vel_all, lambda, ...
        ref_coor0, -10);
    fprintf('[仿真] 杂波生成完毕。\n\n');
end

%% =========================================================================
% PART 6：运动补偿（MOCO）+ 脉冲压缩（参考点 = 蜂群中心 ref_coor0）
% =========================================================================
moco_ref = ref_coor0;                               % 固定场景参考点
dist1 = vecnorm(ant_pos_all  - moco_ref, 2, 2);     % 含姿态扰动瞬时斜距
dist2 = vecnorm(ant_ideal_all - moco_ref, 2, 2);    % 理想瞬时斜距
delta_swing_all     = dist1 - dist2;
swing_phase_2d_corr = exp(1j * 4*pi * delta_swing_all / lambda);

t_ref = (0:N_fast-1) / fs;
s_t   = rectpuls(t_ref - Tc/2, Tc) .* exp(1j*pi*K*t_ref.^2);
X_tx  = fft(s_t, N_fast);

if exist(fullfile(pwd,'snr_override.txt'),'file'), SNR_db = str2double(fileread(fullfile(pwd,'snr_override.txt'))); else, SNR_db = -15; end
s_total      = awgn(s_rx_raw + s_clutter_full, SNR_db);
s_total_moco = s_total .* swing_phase_2d_corr;

s_range_raw  = zeros(Num_chirp, N_fast, 'like', 1+1j);
s_range_moco = zeros(Num_chirp, N_fast, 'like', 1+1j);
for i = 1:Num_chirp
    Y_raw  = fft(s_total(i,:),      N_fast);
    Y_moco = fft(s_total_moco(i,:), N_fast);
    s_range_raw(i,:)  = ifft(Y_raw  .* conj(X_tx) .* win_range);
    s_range_moco(i,:) = ifft(Y_moco .* conj(X_tx) .* win_range);
end
fprintf('[MOCO] 相位补偿 + 脉冲压缩完成。\n\n');

%% =========================================================================
% PART 7：创新点1 — 运动状态驱动的杂波脊线参数（每 CPI 标量预测，供 C 与 IMU 先验）
% =========================================================================
fprintf('[创新1] 运动状态驱动的杂波脊线参数建模...\n');
ant0_approx = pos_all(1,:) + radar_offset;
beam_dir    = atan2(IN_coor0(2)-ant0_approx(2), IN_coor0(1)-ant0_approx(1));
theta_bw    = deg2rad(30);
clutter_model = compute_clutter_ridge_model(vel_all, ant_vel_all, ang_rate_all, ...
    pos_all, ang_all, lambda, PRI, N_slow, step, num_frames, Num_chirp, radar_offset, ...
    beam_dir, theta_bw, IN_coor0, 0);
fprintf('[创新1] 已计算 %d 帧平台几何杂波分量 | σ_c_platform∈[%.2f, %.2f] m/s\n\n', ...
    num_frames, min(clutter_model.sigma_c_est), max(clutter_model.sigma_c_est));

%% =========================================================================
% PART 8：处理配置（融合权重 + 时频处理参数）
% =========================================================================
BETA_FUSION = 0.6;    % IMU 先验权重（RD 软陷波 C 用）
KAPPA_BW    = 2.5;    % 陷波/掩模带宽扩展因子
ALPHA_EMA      = 0.7;       % EMA 帧间平滑（消除瞬时估计噪声，大=慢响应更稳定）
sigma_c_ema    = 1.5;       % 运行时 EMA 状态（初始=保守估计，几帧后收敛到数据）
sigma_c_total_arr = zeros(num_frames,1); % 逐帧总展宽记录（可视化）

% ── 时频处理配置（创新点2/3）──────────────────────────────────────────────
tfc.L         = 64;            % STFT 窗长（慢时间样本）
tfc.H         = 16;            % 帧移（L/4，高重叠保证 WOLA 完美重建）
tfc.Nfft      = N_slow;        % STFT FFT 点数（= N_slow，使 v 轴与 v_tf 一致）
tfc.pad       = tfc.L/2;       % 两端补零，覆盖整个 CPI
tfc.w         = hann(tfc.L);   % 分析=合成窗（Hann）
tfc.v_tf      = v_tf;
tfc.dk_search = 30;            % Viterbi 搜索半窗（bin），≈±6.5 m/s
tfc.sigmaT    = 3;             % 帧间转移平滑标准差（bin）
tfc.lamT      = 0.4;           % 转移项权重
tfc.kappa     = KAPPA_BW;
tfc.Hmin      = 0.01;          % 目标最低通过增益（目标失真约束）
tfc.alpha     = 1.0;           % 陷波中心最大深度
tfc.Bflat_cap = 2.5;           % 平顶半宽上限 (m/s)，防吞目标
DELTA_MARGIN = 0.1;

% 3dB带宽第一次滤波统计跟踪（替换固定比例 β=0.8）
bw_3dB_arr = zeros(num_frames,1);      % 逐帧3dB带宽 (m/s)
bw_3dB_ratio_arr = zeros(num_frames,1);  % bw_3dB / B_flat_final 比值

fprintf('[配置] 融合 β=%.1f, κ=%.1f | STFT: L=%d, H=%d, Nfft=%d, 帧数≈%d\n', ...
    BETA_FUSION, KAPPA_BW, tfc.L, tfc.H, tfc.Nfft, ...
    floor((N_slow+2*tfc.pad-tfc.L)/tfc.H)+1);
fprintf('STFT 窗时长=%.1f ms | 帧间隔=%.1f ms | 窗内速度分辨率≈%.2f m/s\n\n', ...
    tfc.L*PRI*1e3, tfc.H*PRI*1e3, (1/(tfc.L*PRI))*lambda/2);
fprintf('[新增] 第一次滤波改用杂波3dB带宽(数据驱动)，替代固定比例(β=0.8)\n');
fprintf('        δ_margin=%.2f m/s | 最低带宽下限=0.5 m/s\n\n', DELTA_MARGIN);
%% =========================================================================
% PART 9：主处理循环
% =========================================================================
fprintf('[主循环] 开始分帧处理（共 %d 帧）...\n', num_frames);

% 性能记录
SCR_A = zeros(num_frames,1); SCR_B = zeros(num_frames,1);
SCR_C = zeros(num_frames,1); SCR_D = zeros(num_frames,1);
SCR_E = zeros(num_frames,1); SCR_F = zeros(num_frames,1);

TERI_C = zeros(num_frames,1); TERI_D = zeros(num_frames,1);
TERI_E = zeros(num_frames,1); TERI_F = zeros(num_frames,1);
TERI_A = zeros(num_frames,1); TERI_B = zeros(num_frames,1);
tgt_separable = false(num_frames,1);  % true: 所有目标都在杂波带外，TERI 可信

% 逐目标指标存储（用于多目标汇总输出）
SCR_D_tgt = zeros(num_frames, N_tgt);
SCR_D0 = zeros(num_frames,1); TERI_D0 = zeros(num_frames,1);
Pd_slow_D0 = zeros(num_frames,1); Pd_fast_D0 = zeros(num_frames,1); pd_slow_D0_hit = 0; pd_fast_D0_hit = 0;
SCR_Dfx = zeros(num_frames,1); TERI_Dfx = zeros(num_frames,1); pd_slow_Dfx_hit = 0; pd_fast_Dfx_hit = 0;
Pd_slow_Dfx = zeros(num_frames,1); Pd_fast_Dfx = zeros(num_frames,1);
SCR_Dfn = zeros(num_frames,1); TERI_Dfn = zeros(num_frames,1); pd_slow_Dfn_hit = 0; pd_fast_Dfn_hit = 0;
Nn = 5; SCR_Nc = zeros(num_frames, Nn); Pd_slow_Nc = zeros(num_frames, Nn); pd_slow_hit_Nc = zeros(1,Nn); pd_fast_hit_Nc = zeros(1,Nn);
fa_cnt_D0 = 0; fa_cnt_Dfx = 0; fa_cnt_Dfn = 0; Pd_slow_Dfn = zeros(num_frames,1); Pd_fast_Dfn = zeros(num_frames,1);
SCR_A_tgt = zeros(num_frames, N_tgt);
SCR_B_tgt = zeros(num_frames, N_tgt);
SCR_C_tgt = zeros(num_frames, N_tgt);
SCR_E_tgt = zeros(num_frames, N_tgt);
SCR_F_tgt = zeros(num_frames, N_tgt);
TERI_A_tgt = zeros(num_frames, N_tgt);
TERI_B_tgt = zeros(num_frames, N_tgt);
TERI_C_tgt = zeros(num_frames, N_tgt);
TERI_D_tgt = zeros(num_frames, N_tgt);
TERI_E_tgt = zeros(num_frames, N_tgt);
TERI_F_tgt = zeros(num_frames, N_tgt);
tgt_vr_all = zeros(num_frames, N_tgt);   % 逐帧各目标径向速度 (m/s)

pts_A = zeros(num_frames,1); pts_B = zeros(num_frames,1);
pts_C = zeros(num_frames,1); pts_D = zeros(num_frames,1);

BW3_C = zeros(num_frames,1); BW3_D = zeros(num_frames,1);
mu_c_all = zeros(num_frames,1);
sigma_c_all = zeros(num_frames,1);
peak_found = false(num_frames,1);

% 先验引导目标定位参数
W_d_search = 4; W_r_search = 5; gamma_min_dB = 4.0;

% CFAR 参数
Tr_go=8; Gr_go=2; Td_go=4; Gd_go=2; Pfa_work=1e-4;
WIDE_NOTCH_BW = 4.0;   % 方案B 宽陷波半宽 (m/s)
FIX_NOTCH_BW  = 0.5;   % 原方案A 固定零陷半宽 (m/s)，现已替换为 MTI，仅作 MID 调试存档

f_mid = 30;   % 保存用于可视化的中间帧
MID = struct();

% ── ROC 曲线扫参设置 ──
PFA_SWEEP = logspace(-6, -1.5, 25);  % 1e-6 ~ 3.16e-2, 25 pts log-uniform
nPfa = numel(PFA_SWEEP);
N_ref_cfar = 2*(Td_go + Tr_go);      % GO-CFAR: 总参考单元数 = 2*(4+8) = 24
d_cfar_lo = Td_go+Gd_go+1; d_cfar_hi = N_slow-Td_go-Gd_go;
r_cfar_lo = Tr_go+Gr_go+1; r_cfar_hi = N_fast-Tr_go-Gr_go;
n_valid_cfar = (d_cfar_hi-d_cfar_lo+1)*(r_cfar_hi-r_cfar_lo+1);
valid_cfar_mask = false(N_slow, N_fast);
valid_cfar_mask(d_cfar_lo:d_cfar_hi, r_cfar_lo:r_cfar_hi) = true;
pd_count = zeros(nPfa, 6);     % [nPfa × 方案A/B/C/D] 目标命中计数
fa_roc_cnt = zeros(nPfa, 6);   % [nPfa × 方案A/B/C/D] 虚警计数
pd_sep_cnt = zeros(nPfa, 6);   % 目标可分离帧的命中计数
pd_clut_cnt = zeros(nPfa, 6);  % 目标在杂波带内帧的命中计数
n_sep_frames = 0;               % 可分离总帧数（目标在 CFAR 有效区）
n_clut_frames = 0;              % 杂波带内总帧数

% ââ éç®æ æ£æµç»è®¡ï¼æ¢é/å¿«éç®æ ä¿æçï¼ââ
n_slow_tgtframes = 0;               % æ¢éç®æ -å¸§å¯¹æ»æ°ï¼in-clutterï¼
n_fast_tgtframes = 0;               % å¿«éç®æ -å¸§å¯¹æ»æ°ï¼out-of-clutterï¼
pd_slow_cnt = zeros(nPfa, 6);       % [nPfa x 4] æ¢éç®æ æ£æµæ¬¡æ°
pd_fast_cnt = zeros(nPfa, 6);       % [nPfa x 4] å¿«éç®æ æ£æµæ¬¡æ°

for f_idx = 1:num_frames
    s_idx = (f_idx-1)*step + 1;
    e_idx = s_idx + N_slow - 1;
    if e_idx > Num_chirp; break; end
    frame_idx = s_idx:e_idx;

    % ── 帧中心 IMU 姿态/位置 ──
    m_fc = min(max(round((s_idx+e_idx)/2),1), Num_chirp);
    t_fc = (m_fc-1)*PRI;
    R_fc = build_R_zyx_local(deg2rad(ang_all(m_fc,1)), deg2rad(ang_all(m_fc,2)), deg2rad(ang_all(m_fc,3)));
    ant_pos_fc = pos_all(m_fc,:) + (R_fc*radar_offset.').';
    slow_idx = (0:N_slow-1).';

    % ── 平台平动多普勒补偿（IMU）：把杂波中心拉回 ~0 ──
    e_los_bs   = R_fc;
    v_rad_plat = mean(vel_all(frame_idx,:) * e_los_bs, 1);
    fd_shift   = 2 * norm(v_rad_plat) / lambda;
    shift_plat = exp(-1j * 2*pi * fd_shift * slow_idx * PRI);

    frm_raw  = s_range_raw(frame_idx,:);
    frm_moco = s_range_moco(frame_idx,:);
    frm_raw_comp  = frm_raw  .* shift_plat;          % 仅平动补偿（方案A用）
    frm_moco_comp = frm_moco .* shift_plat;          % 平动+姿态MOCO（B/C/D用）

    % ── 共用：MOCO 后未陷波 RD（残余杂波集中在 ~0）──
    rd_pre = fftshift(fft(frm_moco_comp .* win_doppler, N_slow, 1), 1);

    % ── 杂波中心/谱宽融合估计（创新点1+2 的标量部分）──
    mu_c_att_prior = clutter_model.mu_c_att(f_idx);
    % 【断裂点2 修正】先确定杂波主导距离区（近距离杂波区），只在该区估计 σ_c，
    %   避免远距离纯噪声门把 RD 域杂波谱峰稀释、导致 σ_c 偏小。此处的距离窗与
    %   方案 D 后续的 tfc.clut_bins 保持一致（近距离 50 起，远端留出目标余量）。
    clut_lo_est = 50;
    clut_hi_est = min(N_fast, round(0.6*N_fast));   % 近~中距离杂波区，不依赖目标真值
    clut_bins_est = clut_lo_est:clut_hi_est;
    %%%%============================================================
    [mu_c_data, sigma_c_data] = estimate_clutter_stats_rd(rd_pre, v_tf, mu_c_att_prior, 4.0, clut_bins_est); % 仅杂波距离区估计残余杂波中心与谱宽
    %%%%============================================================
    mu_c_imu    = mu_c_att_prior;
    sigma_c_imu = clutter_model.sigma_c_est(f_idx);
    mu_c_final    = BETA_FUSION*mu_c_imu + (1-BETA_FUSION)*mu_c_data;
    mu_c_final    = max(min(mu_c_final, 1.5), -1.5);  %杂波谱中心
    sigma_c_final = max(sigma_c_imu + sigma_c_data, 0.3); %杂波谱展宽
    B_flat_final  = min(KAPPA_BW*sigma_c_final, tfc.Bflat_cap); %tfc.Bflat_cap，平顶半宽上限
    B_roll_final  = 0.5*sigma_c_final;
    mu_c_all(f_idx) = mu_c_final;
    sigma_c_all(f_idx) = sigma_c_final;

    % ── 逐脉冲姿态径向速度（D 方案蜿蜒先验）──
    v_mean_f     = mean(vel_all(frame_idx,:),1);
    e_bore       = v_mean_f/(norm(v_mean_f)).'; %方向向量
    v_arm_f      = ant_vel_all(frame_idx,:) - vel_all(frame_idx,:);
    dv_att_pulse = -(v_arm_f .* e_bore);             % [N_slow×1]

    % ── 数据驱动总展宽估计（核心：用数据代替固定 σ_icm）──
    %  σ_c_platform = sqrt(σ_beam² + σ_att²)：由 IMU 精确给出，σ_icm 未知但不重要；
    %  σ_c_data：从回波实测（自动包含 σ_icm 的实际影响）；
    %  结果：max(先验下界, 数据实测) → 陷波宽度紧贴真实杂波边界，最大保留边缘目标

    sigma_c_platform = clutter_model.sigma_c_est(f_idx);          % 平台运动引起的先验下界
    sigma_c_meas     = max(sigma_c_platform, sigma_c_data);        % 数据主导，先验兜底
    sigma_c_ema      = ALPHA_EMA*sigma_c_ema + (1-ALPHA_EMA)*sigma_c_meas;  % EMA 平滑
    sigma_c_ema      = max(sigma_c_ema, 0.8 * sigma_c_platform);   % 避免低估
    sigma_c_total_arr(f_idx) = sigma_c_ema;
    % 更新带宽参数（供 C/D 共用）
    sigma_c_imu    = sigma_c_ema;
    sigma_c_final  = max(sigma_c_ema, 0.3);
    sigma_c_all(f_idx) = sigma_c_final;

    % ── 多目标真值位置（仅评估用）──
    tgt_bin_arr = zeros(N_tgt,1);
    tgt_vr_arr  = zeros(N_tgt,1);
    tgt_dop_guess_arr = zeros(N_tgt,1);
    tgt_R_arr = zeros(N_tgt,1);
    for k = 1:N_tgt
        tp_k    = tgt(k).coor0 + t_fc*tgt(k).vt;
        los_vec = tp_k - ant_pos_fc;
        e_los_n = (los_vec/norm(los_vec)).';
        tgt_R_arr(k) = norm(los_vec);
        tgt_bin_arr(k) = min(max(round(tgt_R_arr(k)/range_res) + 1, 1), N_fast);
        tgt_vr_arr(k)  = dot(tgt(k).vt, e_los_n.');
        [~, tgt_dop_guess_arr(k)] = min(abs(v_tf - (-tgt_vr_arr(k))));
    end
    % 最近目标作为 clut_hi / tf_ridge_filter / 主参考
    [~, idx_closest] = min(tgt_R_arr);
    tgt_bin = tgt_bin_arr(idx_closest);
    tgt_vr  = tgt_vr_arr(idx_closest);
    tgt_dop_guess = tgt_dop_guess_arr(idx_closest);
    tgt_R_n = tgt_R_arr(idx_closest);

    % ── 方案 A：仅平动补偿 + 传统 MTI 滤波（二脉冲对消）──
    % 传统 MTI 原理：相邻两个脉冲作差分，抑制静止/慢速杂波（响应在零频处为零）。
    % 滤波器系数 b=[1,-1]，沿慢时间维（dim=1）逐距离门滤波；
    % 首脉冲无前驱数据，置零以消除初始瞬态，再施加多普勒窗后做 FFT。
    b_mti = [1, -1];                                      % 二脉冲对消系数
    frm_A_mti = filter(b_mti, 1, frm_raw_comp, [], 1);   % 沿慢时间维 (dim=1) 滤波
    frm_A_mti(1, :) = 0;                                  % 首脉冲无前驱，清零初始瞬态
    rd_A = fftshift(fft(frm_A_mti .* win_doppler, N_slow, 1), 1);

    % ── 方案 B：MOCO + 宽硬陷波（过抑制基线）──
    H_B  = make_wide_notch(v_tf, mu_c_final, WIDE_NOTCH_BW, N_slow, N_fast);
    rd_B = rd_pre .* H_B;

    % ── 方案 C：MOCO + RD 域自适应软陷波（原 v2 方法）──
    H_C  = adaptive_soft_notch(v_tf, mu_c_final, B_flat_final, B_roll_final, ...
        tfc.alpha, tfc.Hmin, N_slow, N_fast);
    rd_C0 = rd_pre .* H_C;

    % ── 候选目标检测 → 2D 局部保护单元（C 与 D 共用）──
    %  【关键修正】早期版本把保护施加到"整条多普勒行"(H_C(dop,:)=1)，会把残余的
    %  零频杂波整行增益抬回 1 → 形成 v≈0 贯穿全距离的亮线。此处改为：
    %    (1) 只保护"孤立点目标"附近的 2D 局部单元，绝不整行；
    %    (2) 排除落在陷波带内(|v-μ_c|<=B_flat)的候选——它们基本是残余杂波而非目标。
    cand_mask = detect_candidate_targets(abs(rd_C0).^2, 12);
    protect_cells = zeros(0,2);                 % [d_idx, r_idx]
    cl = find(cand_mask).';
    for k = cl
        [d_idx, r_idx] = ind2sub([N_slow, N_fast], k);
        in_notch = abs(v_tf(d_idx) - mu_c_final) <= B_flat_final;
        if ~in_notch && isolated_peak_test(abs(rd_C0).^2, r_idx, d_idx, 5, 10)
            protect_cells(end+1,:) = [d_idx, r_idx];   %#ok<SAGROW>
        end
    end

    H_C2 = H_C;                                  % 2D 局部保护
    for ci = 1:size(protect_cells,1)
        dd = protect_cells(ci,1); rr = protect_cells(ci,2);
        H_C2(max(1,dd-2):min(N_slow,dd+2), max(1,rr-3):min(N_fast,rr+3)) = 1.0;
    end
    rd_C = rd_pre .* H_C2;

    % ── 方案 D：MOCO + 慢时间STFT + Viterbi脊线追踪 + 目标保持型时频掩模 ──
    % 杂波主导的距离区间（排除目标附近），用于脊线追踪能量积累
    clut_lo = 50;
    Nc_list = [5 10 20 30 50];   % 距离单元数敏感性扫描（杂波距离区单元数）
    clut_hi = max(clut_lo+1, min(tgt_bin_arr)-30);
    tfc.clut_bins = clut_lo:min(clut_hi, N_fast);

    % 从RD域数据估计杂波峰值的3dB带宽（替代固定比例 β=0.8）
    bw_3dB = compute_clutter_3dB_bw(rd_pre, v_tf, mu_c_final, clut_bins_est);
    B_core_f = max(bw_3dB, 0.5);  % 最低 0.5 m/s 保证基本压制
    bw_3dB_arr(f_idx) = bw_3dB;
    bw_3dB_ratio_arr(f_idx) = bw_3dB / max(B_flat_final, eps);
    rd_pass1 = apply_conservative_notch(rd_pre, v_tf, mu_c_final, B_core_f);

    % 时频掩模本身不做整行保护(掩模对带外目标增益≈1)，输出 RD 后再做 2D 局部恢复
    [rd_D_pure, diagD] = tf_ridge_filter(frm_moco_comp, win_doppler, tfc, ...
        dv_att_pulse, sigma_c_imu, tgt_bin);

    cand_mask_d = detect_candidate_targets(abs(rd_pass1).^2,10);
    protect_cells_D = zeros(0,2);
    n_cand_raw_D  = 0;   % 候选检测原始数（用于诊断）
    n_cand_kept_D = 0;   % 速度门控后保留数

    cl = find(cand_mask_d).';

    B_flat_for_gate = max(B_flat_final, diagD.B_flat_used_ms);
    for k = cl
        [d_idx, r_idx] = ind2sub([N_slow, N_fast], k);
        n_cand_raw_D = n_cand_raw_D + 1;
        v_cand = v_tf(d_idx);

        % 速度域门控：杂波带内（含安全余量）的候选直接拒绝
        in_clutter_zone = abs(v_cand - mu_c_final) <= (B_flat_for_gate + DELTA_MARGIN);

        % 双重验证：速度在杂波带外 AND 孤立点测试通过
        if ~in_clutter_zone && isolated_peak_test(abs(rd_pass1).^2, r_idx, d_idx, 4, 8)
            protect_cells_D(end+1,:) = [d_idx, r_idx];   %#ok<SAGROW>
            n_cand_kept_D = n_cand_kept_D + 1;
        end
    end

    rd_D = rd_D_pure;
    for ci = 1:size(protect_cells_D,1)
        dd = protect_cells_D(ci,1);
        rr = protect_cells_D(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2);
        rrg = max(1,rr-3):min(N_fast,rr+3);
        % 从未陷波的原始谱中恢复目标能量，大幅提高 TERI
        rd_D(dr, rrg) = rd_pre(dr, rrg);
    end

    rd_D0 = rd_D_pure;    % D without Phase-3 (no target-energy compensation)
    [rd_Dfx_pure, ~] = tf_ridge_filter(frm_moco_comp, win_doppler, tfc, dv_att_pulse, sigma_c_imu, tgt_bin, 4.0);  % w/o adaptive mask: FIXED wide (4.0 m/s)
    rd_Dfx = rd_Dfx_pure;
    for ci = 1:size(protect_cells_D,1)
        dd = protect_cells_D(ci,1); rr = protect_cells_D(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2); rrg = max(1,rr-3):min(N_fast,rr+3);
        rd_Dfx(dr, rrg) = rd_pre(dr, rrg);
    end
    [rd_Dfn_pure, ~] = tf_ridge_filter(frm_moco_comp, win_doppler, tfc, dv_att_pulse, sigma_c_imu, tgt_bin, 0.5);  % w/o adaptive mask: FIXED narrow (0.5 m/s)
    rd_Dfn = rd_Dfn_pure;
    for ci = 1:size(protect_cells_D,1)
        dd = protect_cells_D(ci,1); rr = protect_cells_D(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2); rrg = max(1,rr-3):min(N_fast,rr+3);
        rd_Dfn(dr, rrg) = rd_pre(dr, rrg);
    end
    % ── 方案 E：MOCO + 子空间/特征对消（滑动子CPI的 RPCA/eigencanceler）──
    Ysub = frm_moco_comp(:, tfc.clut_bins);                 % 杂波区慢时间快拍 [N_slow×Nc]
    Rc = (Ysub*Ysub')/max(size(Ysub,2),1);                  % 杂波协方差
    [Vc,Dc] = eig((Rc+Rc')/2); ev = real(diag(Dc));
    [ev,ord] = sort(ev,'descend'); Vc = Vc(:,ord);
    Kc = sum(ev > 0.05*ev(1)); Kc = min(max(Kc,1), N_slow-2);
    Uc = Vc(:,1:Kc);                                        % 杂波子空间
    frm_E = (eye(N_slow) - Uc*Uc') * frm_moco_comp;         % 正交补投影 → 对消杂波
    rd_E = fftshift(fft(frm_E .* win_doppler, N_slow, 1), 1);

    % ── 方案 F：无 IMU 先验的数据驱动时频脊线追踪（同类方法，去先验）──
    tfc_F = tfc; tfc_F.dk_search = N_slow;                  % 搜索窗放开 → 无先验约束
    [rd_F_pure, ~] = tf_ridge_filter(frm_moco_comp, win_doppler, tfc_F, ...
        zeros(N_slow,1), 0, tgt_bin);                       % 无IMU蜿蜒/无IMU谱宽
    rd_F = rd_F_pure;
    for ci = 1:size(protect_cells_D,1)
        dd = protect_cells_D(ci,1); rr = protect_cells_D(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2); rrg = max(1,rr-3):min(N_fast,rr+3);
        rd_F(dr, rrg) = rd_pre(dr, rrg);                    % 同D的目标能量恢复
    end
% ── 多目标峰值定位 + 性能指标（逐目标）──
    scr_a_k = zeros(N_tgt,1); scr_b_k = zeros(N_tgt,1);
    scr_c_k = zeros(N_tgt,1); scr_d_k = zeros(N_tgt,1); scr_e_k = zeros(N_tgt,1); scr_f_k = zeros(N_tgt,1); scr_d0_k = zeros(N_tgt,1); scr_dfx_k = zeros(N_tgt,1); scr_dfn_k = zeros(N_tgt,1);
    teri_a_k = zeros(N_tgt,1); teri_b_k = zeros(N_tgt,1);
    teri_c_k = zeros(N_tgt,1); teri_d_k = zeros(N_tgt,1); teri_e_k = zeros(N_tgt,1); teri_f_k = zeros(N_tgt,1); teri_d0_k = zeros(N_tgt,1); teri_dfx_k = zeros(N_tgt,1); teri_dfn_k = zeros(N_tgt,1);
    sep_k = false(N_tgt,1);
    peak_ok_k = false(N_tgt,1);
    d_ref_k_arr = zeros(N_tgt,1);
    r_ref_k_arr = zeros(N_tgt,1);

    tgt_zone = false(N_slow, N_fast);   % 多目标联合检测区
    for k = 1:N_tgt
        [d_ref_k, r_ref_k, ~, ~, ok_k] = find_target_peak_prior(rd_pre, ...
            tgt_dop_guess_arr(k), tgt_bin_arr(k), W_d_search, W_r_search, gamma_min_dB);
        peak_ok_k(k) = ok_k;
        if ~ok_k
            d_ref_k = tgt_dop_guess_arr(k); r_ref_k = tgt_bin_arr(k);
        end
        d_ref_k_arr(k) = d_ref_k; r_ref_k_arr(k) = r_ref_k;

        [scr_a_k(k), P_a] = compute_scr(rd_A, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_b_k(k), P_b] = compute_scr(rd_B, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_c_k(k), P_c] = compute_scr(rd_C, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_d_k(k), P_d] = compute_scr(rd_D, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_e_k(k), P_e] = compute_scr(rd_E, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_f_k(k), P_f] = compute_scr(rd_F, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_d0_k(k), P_d0] = compute_scr(rd_D0, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_dfx_k(k), P_dfx] = compute_scr(rd_Dfx, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        [scr_dfn_k(k), P_dfn] = compute_scr(rd_Dfn, d_ref_k, r_ref_k, N_slow, N_fast, clut_bins_est);
        r_teri = max(1,r_ref_k-1):min(N_fast,r_ref_k+1);
        d_teri = max(1,d_ref_k-1):min(N_slow,d_ref_k+1);
        % ── 原始总功率（含目标+杂波+噪声）──
        P_pre = mean(abs(rd_pre(d_teri,r_teri)).^2, 'all') + eps;

        % ═══════════════════════════════════════════════════════════════
        %  思路三改进：扣除杂波基底后的目标能量保持率
        %  在相同多普勒门(d_teri)上，取杂波距离区(clut_bins_est)中
        %  远离当前目标 ±20 门的位置，用功率中位数估计该多普勒的
        %  杂波基底功率 P_clut，然后从总功率中扣除：
        %    P_tgt_净 = max(P_total - P_clut, eps)
        %    TERI     = P_tgt_净_after / P_tgt_净_before
        % ═══════════════════════════════════════════════════════════════
        guard_r = 20;  % 排除目标附近距离门，防止目标能量泄漏污染杂波估计
        clut_ref = clut_bins_est(abs(clut_bins_est - r_ref_k) > guard_r);
        if length(clut_ref) < 10
            clut_ref = clut_bins_est;   % 保底：如果排除后样本太少则全部使用
        end
        % 各方案在目标多普勒处的杂波基底功率（中位数，对 outlier 鲁棒）
        P_clut_pre = median(abs(rd_pre(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_A   = median(abs(rd_A(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_B   = median(abs(rd_B(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_C   = median(abs(rd_C(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_D   = median(abs(rd_D(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_E   = median(abs(rd_E(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_F   = median(abs(rd_F(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_D0  = median(abs(rd_D0(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_Dfx = median(abs(rd_Dfx(d_teri, clut_ref)).^2, 'all') + eps;
        P_clut_Dfn = median(abs(rd_Dfn(d_teri, clut_ref)).^2, 'all') + eps;

        % 扣除杂波后的目标净能量
        P_tgt_pre = max(P_pre - P_clut_pre, eps);
        P_tgt_A   = max(P_a  - P_clut_A, eps);
        P_tgt_B   = max(P_b  - P_clut_B, eps);
        P_tgt_C   = max(P_c  - P_clut_C, eps);
        P_tgt_D   = max(P_d  - P_clut_D, eps);
        P_tgt_E   = max(P_e  - P_clut_E, eps);
        P_tgt_F   = max(P_f  - P_clut_F, eps);
        P_tgt_D0   = max(P_d0 - P_clut_D0, eps);
        P_tgt_Dfx  = max(P_dfx - P_clut_Dfx, eps);
        P_tgt_Dfn  = max(P_dfn - P_clut_Dfn, eps);

        % 真正的目标能量保持率
        teri_a_k(k) = P_tgt_A / P_tgt_pre;
        teri_b_k(k) = P_tgt_B / P_tgt_pre;
        teri_c_k(k) = P_tgt_C / P_tgt_pre;
        teri_d_k(k) = P_tgt_D / P_tgt_pre;
        teri_e_k(k) = P_tgt_E / P_tgt_pre;
        teri_f_k(k) = P_tgt_F / P_tgt_pre;
        teri_d0_k(k) = P_tgt_D0 / P_tgt_pre;
        teri_dfx_k(k) = P_tgt_Dfx / P_tgt_pre;
        teri_dfn_k(k) = P_tgt_Dfn / P_tgt_pre;
        sep_k(k) = abs(-tgt_vr_arr(k) - mu_c_final) > B_flat_final + DELTA_MARGIN;
        % 联合检测区：所有目标的邻域并集
        tgt_zone(max(1,d_ref_k-2):min(N_slow,d_ref_k+2), ...
            max(1,r_ref_k-2):min(N_fast,r_ref_k+2)) = true;
    end
    tgt_zone = tgt_zone & valid_cfar_mask;

    % 聚合指标（最差 SCR，平均 TERI）
    SCR_A(f_idx) = min(scr_a_k);
    SCR_B(f_idx) = min(scr_b_k);
    SCR_C(f_idx) = min(scr_c_k);
    SCR_D(f_idx) = min(scr_d_k);
    SCR_E(f_idx) = min(scr_e_k);
    SCR_F(f_idx) = min(scr_f_k);
    SCR_D_tgt(f_idx,:) = scr_d_k.';
    SCR_A_tgt(f_idx,:) = scr_a_k.';
    SCR_B_tgt(f_idx,:) = scr_b_k.';
    SCR_C_tgt(f_idx,:) = scr_c_k.';
    SCR_E_tgt(f_idx,:) = scr_e_k.';
    SCR_F_tgt(f_idx,:) = scr_f_k.';            % 逐目标存储
    TERI_A(f_idx) = mean(teri_a_k);
    TERI_B(f_idx) = mean(teri_b_k);
    TERI_C(f_idx) = mean(teri_c_k);
    TERI_D(f_idx) = mean(teri_d_k);
    TERI_E(f_idx) = mean(teri_e_k);
    TERI_F(f_idx) = mean(teri_f_k);
    SCR_D0(f_idx) = min(scr_d0_k); TERI_D0(f_idx) = mean(teri_d0_k);
    SCR_Dfx(f_idx) = min(scr_dfx_k); TERI_Dfx(f_idx) = mean(teri_dfx_k);
    SCR_Dfn(f_idx) = min(scr_dfn_k); TERI_Dfn(f_idx) = mean(teri_dfn_k);
    TERI_D_tgt(f_idx,:) = teri_d_k.';
    TERI_E_tgt(f_idx,:) = teri_e_k.';
    TERI_F_tgt(f_idx,:) = teri_f_k.';          % 逐目标存储
    TERI_C_tgt(f_idx,:) = teri_c_k.';          % 逐目标存储（方案C，用于D/C对比）
    TERI_B_tgt(f_idx,:) = teri_b_k.';          % 逐目标存储（方案B）
    TERI_A_tgt(f_idx,:) = teri_a_k.';          % 逐目标存储（方案A）
    tgt_vr_all(f_idx,:) = tgt_vr_arr.';        % 逐目标径向速度存储（用于实验结果分析
    % 可分离标记：全部可分离才标为 separable
    tgt_separable(f_idx) = all(sep_k);
    % éç®æ å¸§è®¡æ°ï¼sep_k(k)=true ä¸ºå¿«é(å¸¦å¤)ï¼false ä¸ºæ¢é(å¸¦å)
    n_slow_tgtframes = n_slow_tgtframes + sum(~sep_k);
    n_fast_tgtframes = n_fast_tgtframes + sum(sep_k);
    BW3_C(f_idx) = compute_bw3dB(rd_C, tfc.clut_bins, v_tf, N_slow);
    BW3_D(f_idx) = compute_bw3dB(rd_D, tfc.clut_bins, v_tf, N_slow);
    % 总体峰值锁定标记：任一目标锁定即视为有效
    peak_found(f_idx) = any(peak_ok_k);
    % 保留主参考目标（最近目标）用于 MID 兼容
    d_ref = d_ref_k_arr(idx_closest);
    r_ref = r_ref_k_arr(idx_closest);
    ok = peak_ok_k(idx_closest);

    % ── CFAR 检测（向量化，捕获 Z_map 用于 ROC 扫描）──
    % ── 距离单元数 N_c 敏感性：逐 N_c 重算脊线→SCR/Pd ──
    for iNc = 1:Nn
        tfcN = tfc; tfcN.clut_bins = clut_lo:min(clut_lo+Nc_list(iNc)-1, N_fast);
        [rdN, ~] = tf_ridge_filter(frm_moco_comp, win_doppler, tfcN, dv_att_pulse, sigma_c_imu, tgt_bin);
        for ci = 1:size(protect_cells_D,1)
            dd = protect_cells_D(ci,1); rr = protect_cells_D(ci,2);
            dr = max(1,dd-2):min(N_slow,dd+2); rrg = max(1,rr-3):min(N_fast,rr+3);
            rdN(dr, rrg) = rd_pre(dr, rrg);
        end
        scN = zeros(N_tgt,1);
        for kt = 1:N_tgt
            [scN(kt), ~] = compute_scr(rdN, d_ref_k_arr(kt), r_ref_k_arr(kt), N_slow, N_fast, clut_bins_est);
        end
        SCR_Nc(f_idx,iNc) = min(scN);
        [detN, ~] = my_go_cfar(abs(rdN).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
        for kt = 1:N_tgt
            dk = d_ref_k_arr(kt); rk = r_ref_k_arr(kt);
            if dk<1 || rk<1; continue; end
            hitN = any(detN(max(1,dk-2):min(N_slow,dk+2), max(1,rk-2):min(N_fast,rk+2)), 'all');
            if sep_k(kt); pd_fast_hit_Nc(iNc) = pd_fast_hit_Nc(iNc) + hitN; else; pd_slow_hit_Nc(iNc) = pd_slow_hit_Nc(iNc) + hitN; end
        end
        Pd_slow_Nc(f_idx,iNc) = pd_slow_hit_Nc(iNc) / max(n_slow_tgtframes,1);
    end
    [det_A, Z_A] = my_go_cfar(abs(rd_A).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_B, Z_B] = my_go_cfar(abs(rd_B).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_C, Z_C] = my_go_cfar(abs(rd_C).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_D, Z_D] = my_go_cfar(abs(rd_D).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_E, Z_E] = my_go_cfar(abs(rd_E).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_F, Z_F] = my_go_cfar(abs(rd_F).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_D0, ~] = my_go_cfar(abs(rd_D0).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_Dfx, ~] = my_go_cfar(abs(rd_Dfx).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_Dfn, ~] = my_go_cfar(abs(rd_Dfn).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    fa_cnt_D0  = fa_cnt_D0  + sum(det_D0, 'all')  - sum(det_D0(tgt_zone), 'all');
    fa_cnt_Dfx = fa_cnt_Dfx + sum(det_Dfx,'all')  - sum(det_Dfx(tgt_zone), 'all');
    fa_cnt_Dfn = fa_cnt_Dfn + sum(det_Dfn,'all')  - sum(det_Dfn(tgt_zone), 'all');
    pts_A(f_idx) = sum(det_A,'all'); pts_B(f_idx) = sum(det_B,'all');
    pts_C(f_idx) = sum(det_C,'all'); pts_D(f_idx) = sum(det_D,'all');
    pts_E(f_idx) = sum(det_E,'all'); pts_F(f_idx) = sum(det_F,'all');

    % ---- ROC sweep over Pfa: reuse Z_map, vary threshold only ----
    tgt_in_valid = any(tgt_zone, 'all');
    if tgt_in_valid
        if tgt_separable(f_idx)
            n_sep_frames = n_sep_frames + 1;
        else
            n_clut_frames = n_clut_frames + 1;
        end
    end

    Z_list = {Z_A, Z_B, Z_C, Z_D, Z_E, Z_F};
    pwr_list = {abs(rd_A).^2, abs(rd_B).^2, abs(rd_C).^2, abs(rd_D).^2, abs(rd_E).^2, abs(rd_F).^2};
    for s_idx = 1:6
        Z_s = Z_list{s_idx};
        p_s = pwr_list{s_idx};
        for p_idx = 1:nPfa
            alpha_p = N_ref_cfar*(PFA_SWEEP(p_idx)^(-1/N_ref_cfar)-1);
            det_p = (p_s > alpha_p*Z_s) & valid_cfar_mask;
            if tgt_in_valid && any(det_p(tgt_zone), 'all')
                pd_count(p_idx, s_idx) = pd_count(p_idx, s_idx) + 1;
                if tgt_separable(f_idx)
                    pd_sep_cnt(p_idx, s_idx) = pd_sep_cnt(p_idx, s_idx) + 1;
                else
                    pd_clut_cnt(p_idx, s_idx) = pd_clut_cnt(p_idx, s_idx) + 1;
                end
            end
            fa_roc_cnt(p_idx, s_idx) = fa_roc_cnt(p_idx, s_idx) + ...
                sum(det_p,'all') - sum(det_p(tgt_zone),'all');

            % Per-target detection: separate slow (in-clutter) vs fast (out-of-clutter)
            for k_tgt = 1:N_tgt
                dk = d_ref_k_arr(k_tgt); rk = r_ref_k_arr(k_tgt);
                if dk < 1 || rk < 1 || dk > N_slow || rk > N_fast; continue; end
                d_lo = max(1, dk-2); d_hi = min(N_slow, dk+2);
                r_lo = max(1, rk-2); r_hi = min(N_fast, rk+2);
                tgt_hit = any(det_p(d_lo:d_hi, r_lo:r_hi), 'all');
                if sep_k(k_tgt)
                    pd_fast_cnt(p_idx, s_idx) = pd_fast_cnt(p_idx, s_idx) + tgt_hit;
                else
                    pd_slow_cnt(p_idx, s_idx) = pd_slow_cnt(p_idx, s_idx) + tgt_hit;
                end
            end
        end
    end

    % D without Phase-3: per-target detection at the reference operating Pfa
    for k_tgt = 1:N_tgt
        dk = d_ref_k_arr(k_tgt); rk = r_ref_k_arr(k_tgt);
        if dk<1 || rk<1 || dk>N_slow || rk>N_fast; continue; end
        d_lo = max(1,dk-2); d_hi = min(N_slow,dk+2);
        r_lo = max(1,rk-2); r_hi = min(N_fast,rk+2);
        hit_D0 = any(det_D0(d_lo:d_hi, r_lo:r_hi), 'all');
        if sep_k(k_tgt); pd_fast_D0_hit = pd_fast_D0_hit + hit_D0; else; pd_slow_D0_hit = pd_slow_D0_hit + hit_D0; end
    end
    Pd_slow_D0(f_idx) = pd_slow_D0_hit / max(n_slow_tgtframes,1);
    Pd_fast_D0(f_idx) = pd_fast_D0_hit / max(n_fast_tgtframes,1);
    % D without adaptive mask: per-target detection at the reference operating Pfa
    for k_tgt = 1:N_tgt
        dk = d_ref_k_arr(k_tgt); rk = r_ref_k_arr(k_tgt);
        if dk<1 || rk<1 || dk>N_slow || rk>N_fast; continue; end
        d_lo = max(1,dk-2); d_hi = min(N_slow,dk+2);
        r_lo = max(1,rk-2); r_hi = min(N_fast,rk+2);
        hit_Dfx = any(det_Dfx(d_lo:d_hi, r_lo:r_hi), 'all');
        if sep_k(k_tgt); pd_fast_Dfx_hit = pd_fast_Dfx_hit + hit_Dfx; else; pd_slow_Dfx_hit = pd_slow_Dfx_hit + hit_Dfx; end
    end
    Pd_slow_Dfx(f_idx) = pd_slow_Dfx_hit / max(n_slow_tgtframes,1);
    Pd_fast_Dfx(f_idx) = pd_fast_Dfx_hit / max(n_fast_tgtframes,1);
    % D with fixed narrow mask: per-target detection at the reference operating Pfa
    for k_tgt = 1:N_tgt
        dk = d_ref_k_arr(k_tgt); rk = r_ref_k_arr(k_tgt);
        if dk<1 || rk<1 || dk>N_slow || rk>N_fast; continue; end
        d_lo = max(1,dk-2); d_hi = min(N_slow,dk+2);
        r_lo = max(1,rk-2); r_hi = min(N_fast,rk+2);
        hit_Dfn = any(det_Dfn(d_lo:d_hi, r_lo:r_hi), 'all');
        if sep_k(k_tgt); pd_fast_Dfn_hit = pd_fast_Dfn_hit + hit_Dfn; else; pd_slow_Dfn_hit = pd_slow_Dfn_hit + hit_Dfn; end
    end
    Pd_slow_Dfn(f_idx) = pd_slow_Dfn_hit / max(n_slow_tgtframes,1);
    Pd_fast_Dfn(f_idx) = pd_fast_Dfn_hit / max(n_fast_tgtframes,1);
    % ── 保存中间帧用于可视化（扩展为多目标）──
    if f_idx == f_mid
        MID.rd_A=rd_A; MID.rd_B=rd_B; MID.rd_C=rd_C; MID.rd_D=rd_D; MID.rd_E=rd_E; MID.rd_F=rd_F; MID.rd_pre=rd_pre;
        MID.rd_pass1   = rd_pass1;
        MID.rd_D_pure  = rd_D_pure;
        MID.B_core     = B_core_f;
        MID.bw_3dB     = bw_3dB;
        MID.n_cand_raw = n_cand_raw_D;
        MID.n_cand_kept= n_cand_kept_D;
        MID.diagD=diagD; MID.d_ref=d_ref; MID.r_ref=r_ref; MID.ok=ok;
        MID.mu_c=mu_c_final; MID.B_flat=B_flat_final; MID.B_roll=B_roll_final;
        MID.B_flat_used = diagD.B_flat_used_ms;
        MID.B_flat_gate = B_flat_for_gate;
        MID.protect_cells  = protect_cells;
        MID.protect_cells_D= protect_cells_D;
        MID.protect_cells=protect_cells; MID.tgt_bin=tgt_bin; MID.t_fc=t_fc;
        MID.v_tgt=-tgt_vr; MID.tgt_R=tgt_R_n; MID.sigma_c_ema=sigma_c_ema;
        MID.clut_bins=tfc.clut_bins; MID.WIDE=WIDE_NOTCH_BW; MID.FIX=FIX_NOTCH_BW;
        % ── 多目标扩展信息 ──
        MID.tgt_bin_arr = tgt_bin_arr;
        MID.tgt_vr_arr = tgt_vr_arr;
        MID.tgt_R_arr = tgt_R_arr;
        MID.d_ref_k_arr = d_ref_k_arr;
        MID.r_ref_k_arr = r_ref_k_arr;
        MID.N_tgt = N_tgt;
        MID.scr_per_tgt = struct('A', scr_a_k, 'B', scr_b_k, 'C', scr_c_k, 'D', scr_d_k);
        MID.teri_per_tgt = struct('A', teri_a_k, 'B', teri_b_k, 'C', teri_c_k, 'D', teri_d_k);
        MID.sep_per_tgt = sep_k;
    end

    if mod(f_idx,20)==0
        fprintf('  帧 %3d/%d | μ_c=%.2f | B_3dB=%.2f | B_flat=%.2f | SCR_D_min=%.1f avg=%.1f dB | Tgt_ok=%d/%d\n', ...
            f_idx, num_frames, mu_c_final, bw_3dB, B_flat_final, ...
            min(scr_d_k), mean(scr_d_k), ...
            sum(peak_ok_k), N_tgt);
    end
end
fprintf('\n[主循环] 处理完成。\n\n');

%% =========================================================================
% PART 10：性能指标汇总
% =========================================================================
% 修正 FAR：排除目标邻域 + CFAR 边缘截断校正
[~, pfa_ref_idx] = min(abs(PFA_SWEEP - Pfa_work));
Pd = pd_count(pfa_ref_idx, :) / num_frames;
FAR_A = fa_roc_cnt(pfa_ref_idx,1) / (num_frames * n_valid_cfar);
FAR_B = fa_roc_cnt(pfa_ref_idx,2) / (num_frames * n_valid_cfar);
FAR_C = fa_roc_cnt(pfa_ref_idx,3) / (num_frames * n_valid_cfar);
FAR_D = fa_roc_cnt(pfa_ref_idx,4) / (num_frames * n_valid_cfar);
FAR_E = fa_roc_cnt(pfa_ref_idx,5) / (num_frames * n_valid_cfar);
FAR_F = fa_roc_cnt(pfa_ref_idx,6) / (num_frames * n_valid_cfar);
FAR_D0  = fa_cnt_D0  / (num_frames * n_valid_cfar);   % ablation arms
FAR_Dfx = fa_cnt_Dfx / (num_frames * n_valid_cfar);
FAR_Dfn = fa_cnt_Dfn / (num_frames * n_valid_cfar);
% ROC 曲线数据（全 Pfa 扫描）
Pd_roc  = pd_count / num_frames;
FAR_roc = fa_roc_cnt / (num_frames * n_valid_cfar);
% 分离性 Pd：目标可分离帧 vs 杂波带内帧
n_pd_total_sep = n_sep_frames + n_clut_frames;      % 有效参与 Pd 统计的总帧数
if n_sep_frames > 0
    Pd_sep_roc = pd_sep_cnt / n_sep_frames;
else
    Pd_sep_roc = zeros(nPfa, 6);
end
if n_clut_frames > 0
    Pd_clut_roc = pd_clut_cnt / n_clut_frames;
else
    Pd_clut_roc = zeros(nPfa, 6);
end
% AUC：标准化 ROC 面积（log10(FAR) 积分）
AUC = zeros(1, 6);
for s_idx = 1:6
    x = FAR_roc(:, s_idx);
    y = Pd_roc(:, s_idx);
    vld = isfinite(x) & isfinite(y) & x > 0;
    if sum(vld) >= 2
        AUC(s_idx) = trapz(log10(x(vld)), y(vld)) / ...
            (log10(max(x(vld))) - log10(min(x(vld))));
    end
end
vmask = peak_found;  nv = sum(vmask);

fprintf('=============================================================\n');
fprintf('         性能指标汇总（%d 帧，有效定位帧 %d）\n', num_frames, nv);
fprintf('=============================================================\n');
fprintf('%-26s %9s %9s %9s %9s %9s %9s\n', '指标','A-MTI滤波','B宽陷波','C-RD软陷波','D-时频脊线','E-子空间','F-无IMU脊线');
fprintf('%s\n', repmat('-',1,84));
fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','平均SCR (dB)', ...
    mean(SCR_A), mean(SCR_B), mean(SCR_C), mean(SCR_D), mean(SCR_E), mean(SCR_F));
if nv>0
    fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','平均SCR-有效帧(dB)', ...
        mean(SCR_A(vmask)), mean(SCR_B(vmask)), mean(SCR_C(vmask)), mean(SCR_D(vmask)), mean(SCR_E(vmask)), mean(SCR_F(vmask)));
    fprintf('%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', 'TERI(all frames)', ...
        mean(TERI_A(vmask)), mean(TERI_B(vmask)), mean(TERI_C(vmask)), mean(TERI_D(vmask)), mean(TERI_E(vmask)), mean(TERI_F(vmask)));
    sep_mask = tgt_separable & vmask;
    nv_sep = sum(sep_mask);
    if nv_sep > 5
        fprintf('%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', 'TERI(clean)', ...
            mean(TERI_A(sep_mask)), mean(TERI_B(sep_mask)), ...
            mean(TERI_C(sep_mask)), mean(TERI_D(sep_mask)), mean(TERI_E(sep_mask)), mean(TERI_F(sep_mask)));
    end
end
fprintf('%-26s %9s %9s %9.2f %9.2f %9s %9s\n','杂波-3dB谱宽(m/s)','-','-', ...
    mean(BW3_C(BW3_C>0)), mean(BW3_D(BW3_D>0)), '-', '-');
fprintf('%-26s %9s %9.2f %9.2f %9.2f %9.2f %9.2f\n','SCR增益 vs A (dB)','-', ...
    mean(SCR_B)-mean(SCR_A), mean(SCR_C)-mean(SCR_A), mean(SCR_D)-mean(SCR_A), mean(SCR_E)-mean(SCR_A), mean(SCR_F)-mean(SCR_A));
fprintf('%-26s %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n','修正虚警率 FAR', FAR_A,FAR_B,FAR_C,FAR_D,FAR_E,FAR_F);
fprintf('  [ablation FAR] D=%.3e  D0(noP3)=%.3e  Dfx(fixW)=%.3e  Dfn(fixN)=%.3e\n', FAR_D, FAR_D0, FAR_Dfx, FAR_Dfn);
fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','检测概率 Pd', Pd(1),Pd(2),Pd(3),Pd(4),Pd(5),Pd(6));
if n_pd_total_sep > 0
    pds = Pd_sep_roc(pfa_ref_idx, :);
    pdc = Pd_clut_roc(pfa_ref_idx, :);
    fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','  Pd(带外帧)', pds);
    fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','  Pd(杂波带内)', pdc);
    fprintf('  [可分性Pd] 带外帧=%d | 带内帧=%d | 合计=%d/总帧=%d\n', ...
        n_sep_frames, n_clut_frames, n_pd_total_sep, num_frames);
end
% Per-target slow/fast detection rate (Pfa=1e-4 reference)
Pd_slow = pd_slow_cnt(pfa_ref_idx, :) / max(n_slow_tgtframes, 1);
Pd_fast = pd_fast_cnt(pfa_ref_idx, :) / max(n_fast_tgtframes, 1);
fid_sw = fopen('snr_sweep.csv','a'); fprintf(fid_sw, '%g,%.3f,%.3f,%.3f\n', SNR_db, mean(SCR_D), Pd_slow(4), Pd_fast(4)); fclose(fid_sw);
fprintf('%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','Per-tgt Pd(slow,in-clutter)', Pd_slow);
fprintf('%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','Per-tgt Pd(fast,out-clutter)', Pd_fast);
fprintf('  D vs C slow Pd improvement: %+.1f%%\n', (Pd_slow(4)-Pd_slow(3))*100);
fprintf('  D vs A slow Pd improvement: %+.1f%%\n', (Pd_slow(4)-Pd_slow(1))*100);

fprintf('%-26s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n','AUC (标准化ROC)', AUC);
fprintf('%s\n', repmat('=',1,66));
fprintf('%-26s %9.3f %9s %9s %9s\n','3dB带宽(第一次滤波,m/s)', mean(bw_3dB_arr), '-','-','-');
fprintf('%-26s %9.3f %9s %9s %9s\n','  3dB带宽/B_flat比值', mean(bw_3dB_ratio_arr), '-','-','-');
fprintf('\n[关键对比] 时频脊线法(D) vs RD软陷波(C):\n');
fprintf('  SCR 差异:        %+.2f dB\n', mean(SCR_D)-mean(SCR_C));
if nv>0
    fprintf('  TERI 比值(C+D):  %.3f (D/C, >1 表示D保持更好)\n', ...
        mean(TERI_D(vmask))/(mean(TERI_C(vmask))+eps));
    sep_mask = tgt_separable & vmask;
    if sum(sep_mask) > 5
        fprintf('  TERI(带外,C+D):  %.3f (D/C, 无杂波污染)\n', ...
            mean(TERI_D(sep_mask))/(mean(TERI_C(sep_mask))+eps));
    end
    fprintf('  TERI_A(MTI)=%.3f | TERI_B(宽陷波)=%.3f (基线损伤)\n', ...
        mean(TERI_A(vmask)), mean(TERI_B(vmask)));
end
fprintf('  虚警率差异:      %+.2e\n', FAR_D-FAR_C);
fprintf('  Pd 差异:          %+.2f\n\n', Pd(4)-Pd(3));

%% =========================================================================
% PART 10b：Multi-target summary (swarm detection statistics)
% =========================================================================
fprintf('\n=============================================================\n');
fprintf('         Multi-target Swarm Detection Summary (%d targets)\n', N_tgt);
fprintf('=============================================================\n');
fprintf('%-26s %9s %9s %9s %9s %9s %9s\n', 'Metric','A-MTI','B-WideN','C-RDSoft','D-TFRidge','E-Subspace','F-NoIMU');
fprintf('%s\n', repmat('-',1,84));
fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','Mean SCR_min (dB)', ...
    mean(SCR_A), mean(SCR_B), mean(SCR_C), mean(SCR_D), mean(SCR_E), mean(SCR_F));
fprintf('%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','Mean TERI', ...
    mean(TERI_A), mean(TERI_B), mean(TERI_C), mean(TERI_D), mean(TERI_E), mean(TERI_F));
fprintf('%s\n\n', repmat('=',1,66));

% Per-target SCR stats (Scheme D, valid frames only)
if nv > 5
    fprintf('Per-target SCR_D stats (valid frames only):\n');
    fprintf('  Tgt# | InitPos(x,y,z)    | Vel(vx,vy,vz)      | Speed  | InitRadVel | MeanSCR_D\n');
    fprintf('%s\n', repmat('-',1,90));
    % Compute initial radial velocity for each target
    vr_init = zeros(N_tgt,1);
    ant0 = ant_pos_all(1,:);
    for k = 1:N_tgt
        los0 = tgt(k).coor0 - ant0;
        e_los0 = los0 / norm(los0);
        vr_init(k) = dot(tgt(k).vt, e_los0);
    end
    for k = 1:N_tgt
        scr_k_mean = mean(SCR_D_tgt(vmask, k), 1);
        tgt_speed = norm(tgt(k).vt);
        fprintf('  #%d | (%4d,%4d,%4d) | (%5.1f,%5.1f,%4.1f) | %5.2f | %+8.2f | %7.1f dB\n', ...
            k, round(tgt(k).coor0), tgt(k).vt, ...
            tgt_speed, vr_init(k), scr_k_mean);
    end
    fprintf('\n');
end

%% =========================================================================
% PART 10c：逐目标能量保持率(TERI)对比分析（全方案 A/B/C/D）
% =========================================================================
%  对每个目标，计算各方案滤波后3x3邻域能量/原始RD能量的比率(TERI)，
%  取所有有效帧的平均值，得到每个目标在各方案下的能量保持率。
fprintf('\n=============================================================\n');
fprintf('  逐目标能量保持率(TERI)对比 — 各方案滤波后目标能量/原始RD能量\n');
fprintf('=============================================================\n');

teri_mean_A = zeros(N_tgt, 1);
teri_mean_B = zeros(N_tgt, 1);
teri_mean_C = zeros(N_tgt, 1);
teri_mean_D = zeros(N_tgt, 1);

for k = 1:N_tgt
    teri_mean_A(k) = mean(TERI_A_tgt(vmask, k));
    teri_mean_B(k) = mean(TERI_B_tgt(vmask, k));
    teri_mean_C(k) = mean(TERI_C_tgt(vmask, k));
    teri_mean_D(k) = mean(TERI_D_tgt(vmask, k));
end

% ── 新增：慢/快目标分组的 SCR 与 TERI（六方案）──
scr_mean = zeros(6, N_tgt); teri_mean_all = zeros(6, N_tgt);
scr_mean(1,:) = mean(SCR_A_tgt(vmask,:),1); scr_mean(2,:) = mean(SCR_B_tgt(vmask,:),1);
scr_mean(3,:) = mean(SCR_C_tgt(vmask,:),1); scr_mean(4,:) = mean(SCR_D_tgt(vmask,:),1);
scr_mean(5,:) = mean(SCR_E_tgt(vmask,:),1); scr_mean(6,:) = mean(SCR_F_tgt(vmask,:),1);
teri_mean_all(1,:) = mean(TERI_A_tgt(vmask,:),1); teri_mean_all(2,:) = mean(TERI_B_tgt(vmask,:),1);
teri_mean_all(3,:) = mean(TERI_C_tgt(vmask,:),1); teri_mean_all(4,:) = mean(TERI_D_tgt(vmask,:),1);
teri_mean_all(5,:) = mean(TERI_E_tgt(vmask,:),1); teri_mean_all(6,:) = mean(TERI_F_tgt(vmask,:),1);
teri_mean_E = teri_mean_all(5,:).'; teri_mean_F = teri_mean_all(6,:).';
SL = 1:5; FA = 6:N_tgt;
fprintf('\n=============================================================\n');
fprintf('  慢/快目标分组汇总（慢=T1~T5, 快=T6~T10；六方案）\n');
fprintf('=============================================================\n');
fprintf('%-24s %9s %9s %9s %9s %9s %9s\n','指标','A-MTI','B-WideN','C-RDSoft','D-TFRidge','E-Subspace','F-NoIMU');
fprintf('%-24s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','SCR 慢目标均值(dB)', mean(scr_mean(:,SL),2)');
fprintf('%-24s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','SCR 快目标均值(dB)', mean(scr_mean(:,FA),2)');
fprintf('%-24s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','TERI 慢目标均值', mean(teri_mean_all(:,SL),2)');
fprintf('%-24s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','TERI 快目标均值', mean(teri_mean_all(:,FA),2)');
M_gf = [mean(scr_mean(:,SL),2)'; mean(scr_mean(:,FA),2)'; mean(teri_mean_all(:,SL),2)'; mean(teri_mean_all(:,FA),2)'];
try; writematrix(M_gf, 'EF_slowfast_grouped.csv'); catch; end
fprintf('  Tgt# | Speed |类别  |  A-MTI  | B-WideN | C-RDSoft| D-TFRidge| E-Subsp | F-NoIMU\n');
fprintf('%s\n', repmat('-', 1, 94));
for k = 1:N_tgt
    if k <= 5
        cat = '慢速';
    elseif k <= 7
        cat = '中速';
    else
        cat = '快速';
    end
    tgt_speed = norm(tgt(k).vt);
    fprintf('  #%2d  | %4.2f | %s | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f\n', ...
        k, tgt_speed, cat, ...
        teri_mean_A(k), teri_mean_B(k), teri_mean_C(k), teri_mean_D(k), teri_mean_E(k), teri_mean_F(k));
end
fprintf('%s\n', repmat('-', 1, 94));
fprintf('  平均  |       | 全  | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f\n', ...
    mean(teri_mean_A), mean(teri_mean_B), ...
    mean(teri_mean_C), mean(teri_mean_D), mean(teri_mean_E), mean(teri_mean_F));
fprintf('  慢均  |       | T1~5| %7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f\n', ...
    mean(teri_mean_A(1:5)), mean(teri_mean_B(1:5)), ...
    mean(teri_mean_C(1:5)), mean(teri_mean_D(1:5)), mean(teri_mean_E(1:5)), mean(teri_mean_F(1:5)));
fprintf('  快均  |       |T6~10| %7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f\n', ...
    mean(teri_mean_A(6:10)), mean(teri_mean_B(6:10)), ...
    mean(teri_mean_C(6:10)), mean(teri_mean_D(6:10)), mean(teri_mean_E(6:10)), mean(teri_mean_F(6:10)));
fprintf('%s\n\n', repmat('=', 1, 94));

%% =========================================================================
% PART 10d：慢速目标(T1~T5)保持率分析
% =========================================================================
%  从两个维度评估各方案对慢速目标的保持能力：
%    (a) 能量保持率 — 各方案在T1~T5上的平均TERI
%    (b) 目标保持计数 — TERI > 阈值(0.5)的目标数，反映该方案能"保住"几个慢速目标
%    (c) 慢速目标各方案间的TERI改善倍数
fprintf('=============================================================\n');
fprintf('  慢速目标(T1~T5)保持率分析 — 各方案对5个慢速目标的保持能力\n');
fprintf('=============================================================\n');

SLOW_TGT_IDS = 1:5;
N_slow_tgt = length(SLOW_TGT_IDS);
TERI_KEEP_THRESH = 0.5;   % TERI > 0.5 认为目标能量被有效保持

% (a) 慢速目标能量保持率详情
fprintf('\n(a) 慢速目标(T1~T5)平均能量保持率(TERI)详情：\n');
fprintf('  %-12s %9s %9s %9s %9s %9s %9s\n', '目标', 'A-MTI', 'B-WideN', 'C-RDSoft', 'D-TFRidge', 'E-Subspace', 'F-NoIMU');
fprintf('%s\n', repmat('-', 1, 79));
for ki = 1:N_slow_tgt
    k = SLOW_TGT_IDS(ki);
    fprintf('  T#%-3d(spd=%.2f) %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', ...
        k, norm(tgt(k).vt), ...
        teri_mean_A(k), teri_mean_B(k), teri_mean_C(k), teri_mean_D(k), teri_mean_E(k), teri_mean_F(k));
end
fprintf('%s\n', repmat('-', 1, 79));
fprintf('  慢速平均         %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', ...
    mean(teri_mean_A(SLOW_TGT_IDS)), mean(teri_mean_B(SLOW_TGT_IDS)), ...
    mean(teri_mean_C(SLOW_TGT_IDS)), mean(teri_mean_D(SLOW_TGT_IDS)), mean(teri_mean_E(SLOW_TGT_IDS)), mean(teri_mean_F(SLOW_TGT_IDS)));
fprintf('  快速(T6~10)平均  %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', ...
    mean(teri_mean_A(6:10)), mean(teri_mean_B(6:10)), ...
    mean(teri_mean_C(6:10)), mean(teri_mean_D(6:10)), mean(teri_mean_E(6:10)), mean(teri_mean_F(6:10)));

% (b) 慢速目标保持计数
keep_count_A = sum(teri_mean_A(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
keep_count_B = sum(teri_mean_B(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
keep_count_C = sum(teri_mean_C(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
keep_count_D = sum(teri_mean_D(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
keep_count_E = sum(teri_mean_E(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
keep_count_F = sum(teri_mean_F(SLOW_TGT_IDS) > TERI_KEEP_THRESH);

fprintf('\n(b) 慢速目标保持计数（TERI > %.1f 的目标数 / 共%d个）：\n', ...
    TERI_KEEP_THRESH, N_slow_tgt);
fprintf('  %-12s %9s %9s %9s %9s %9s %9s\n', '方案', 'A-MTI', 'B-WideN', 'C-RDSoft', 'D-TFRidge', 'E-Subspace', 'F-NoIMU');
fprintf('%s\n', repmat('-', 1, 79));
fprintf('  保持数/总  %6d/%-2d %6d/%-2d %6d/%-2d %6d/%-2d %6d/%-2d %6d/%-2d\n', ...
    keep_count_A, N_slow_tgt, keep_count_B, N_slow_tgt, ...
    keep_count_C, N_slow_tgt, keep_count_D, N_slow_tgt, keep_count_E, N_slow_tgt, keep_count_F, N_slow_tgt);
fprintf('  保持率      %8.1f%% %8.1f%% %8.1f%% %8.1f%% %8.1f%% %8.1f%%\n', ...
    keep_count_A/N_slow_tgt*100, keep_count_B/N_slow_tgt*100, ...
    keep_count_C/N_slow_tgt*100, keep_count_D/N_slow_tgt*100, keep_count_E/N_slow_tgt*100, keep_count_F/N_slow_tgt*100);

% 被保持的目标编号
kept_A = find(teri_mean_A(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
kept_B = find(teri_mean_B(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
kept_C = find(teri_mean_C(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
kept_D = find(teri_mean_D(SLOW_TGT_IDS) > TERI_KEEP_THRESH);
fprintf('  被保持目标: A={%s} | B={%s} | C={%s} | D={%s}\n', ...
    sprintf('%d ', kept_A), sprintf('%d ', kept_B), ...
    sprintf('%d ', kept_C), sprintf('%d ', kept_D));

% (c) 慢速目标TERI改善倍数（D vs C, D vs A）
fprintf('\n(c) 慢速目标TERI改善倍数（方案D相对于其他方案）：\n');
fprintf('  %-12s %9s %9s\n', '目标', 'D/C改善', 'D/A改善');
fprintf('%s\n', repmat('-', 1, 34));
for ki = 1:N_slow_tgt
    k = SLOW_TGT_IDS(ki);
    improv_DC = teri_mean_D(k) / max(teri_mean_C(k), eps);
    improv_DA = teri_mean_D(k) / max(teri_mean_A(k), eps);
    fprintf('  T#%-3d(spd=%.2f) %8.2fx %8.2fx\n', ...
        k, norm(tgt(k).vt), improv_DC, improv_DA);
end
fprintf('%s\n', repmat('-', 1, 34));
fprintf('  平均改善       %8.2fx %8.2fx\n', ...
    mean(teri_mean_D(SLOW_TGT_IDS)) / max(mean(teri_mean_C(SLOW_TGT_IDS)), eps), ...
    mean(teri_mean_D(SLOW_TGT_IDS)) / max(mean(teri_mean_A(SLOW_TGT_IDS)), eps));
fprintf('%s\n\n', repmat('=', 1, 55));

%% =========================================================================
% PART 11：可视化
% =========================================================================
t_frames = (1:num_frames)*step*PRI;

% ── 图1：时频面脊线追踪（创新点1/2 核心可视化）──────────────────────────────
if isfield(MID,'diagD')
    dg = MID.diagD;
    Ptf_dB = 10*log10(dg.Ptf/(max(dg.Ptf(:))+eps)+eps);
    tf_t = ((1:dg.Nf)-1)*tfc.H*PRI*1e3;   % STFT 帧时间 (ms, CPI 内)
   
    figure(1);
    imagesc(tf_t, v_tf, Ptf_dB); axis xy; colorbar; clim([-40,0]);
    colormap(gca,parula); hold on;
    plot(tf_t, v_tf(dg.k_pred),  'w--','LineWidth',1.6,'DisplayName','IMU Prior Spinal Line');
    plot(tf_t, dg.v_ridge,       'r-', 'LineWidth',2.0,'DisplayName','Viterbi Ridge Tracking');
    plot(tf_t, dg.v_ridge+dg.B_flat_v, 'm:','LineWidth',1.2,'DisplayName','Flat-top Heavens');
    plot(tf_t, dg.v_ridge-dg.B_flat_v, 'm:','LineWidth',1.2);
    yline(MID.v_tgt, 'g-', sprintf('Desired speed %.1f m/s ', MID.v_tgt, MID.tgt_R), ...
        'LineWidth',1.4,'LabelHorizontalAlignment','left','FontSize',8);
    xlabel('CPI Internal Time (ms)', 'FontSize', 16);
    ylabel('Speed (m/s)', 'FontSize', 16);
    title(sprintf('Clutter time-frequency spectrum and ridge tracking (Frame %d)', f_mid));
    legend('Location','best'); ylim([min(v_tf), max(v_tf)]); grid on;
  
    figure(2);
   
    imagesc(tf_t, v_tf, dg.G); axis xy; colorbar; clim([0,1]);
    colormap(gca, flipud(gray)); hold on;
    plot(tf_t, dg.v_ridge,'r-','LineWidth',1.5,'DisplayName','Spinal line');
    xlabel('CPI Internal Time', 'FontSize', 16);
    ylabel('Speed (m/s)', 'FontSize', 16);
    title('Target-preserving time-frequency mask G(t,v)');
    legend('Location','best'); grid on;
%     sgtitle('创新点1/2：慢时间时频动态杂波建模与 IMU 约束 Viterbi 脊线追踪', ...
%         'FontSize',12,'FontWeight','bold');
end

%% ── 图2a~2e：原始MOCO + 四方案RD 分别展示 ────────────────────────────────
if isfield(MID,'rd_A')
    clim_rd = [-10, 100];
    rd_list = {MID.rd_pre, MID.rd_A, MID.rd_B, MID.rd_C, MID.rd_D, MID.rd_E, MID.rd_F};
    rd_names = {'Raw RD (after MOCO, unfiltered)', ...
        'A：Traditional MTI Filtering', ...
        sprintf('B：MOCO+Wide Notch(±%.1fm/s)',MID.WIDE), ...
        'C：MOCO+RD Adaptive Soft Notch', ...
        'D：MOCO+Time-Frequency Ridge Tracking', ...
        'E：Subspace Clutter Cancellation (RPCA/eigencanceler)', ...
        'F：Data-driven TF Ridge Tracking (no IMU prior)'};
    fig_labels = {'fig2a_rd_orig', 'fig2b_rd_mti', 'fig2c_rd_widenotch', ...
        'fig2d_rd_softnotch', 'fig2e_rd_tfridge', 'fig2f_rd_subspace', 'fig2g_rd_noimu'};

    for s = 1:7
        figure('Name',rd_names{s});
        imagesc(range_axis, v_tf, 20*log10(abs(rd_list{s})+eps)); axis xy; colorbar;
        clim(clim_rd); hold on;
        % 标注所有目标位置（红色圆圈）
        if isfield(MID,'d_ref_k_arr') && isfield(MID,'r_ref_k_arr')
            colors = {'ro','ro','ro','ro','ro','ro','ro','ro','ro','ro'};
            for kk = 1:min(MID.N_tgt, length(colors))
                plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)), ...
                    colors{kk},'MarkerSize',8,'LineWidth',1.5);
            end
        else
            plot(range_axis(MID.r_ref), v_tf(MID.d_ref), 'r^','MarkerSize',12,'LineWidth',2);
        end
        % 标注杂波中心线 μ_c（原始RD和方案D标注）
        if s==1 || s==6
            yline(MID.mu_c, 'm:', '\mu_c','LineWidth',1.0);
        end
        % 方案D额外标注保护单元
        if s==6
            if isfield(MID,'protect_cells_D') && ~isempty(MID.protect_cells_D)
                plot(range_axis(MID.protect_cells_D(:,2)), v_tf(MID.protect_cells_D(:,1)), ...
                    'gs','MarkerSize',7,'LineWidth',1.2,'DisplayName','Protection Unit');
            end
        end
        xlabel('range (m)'); ylabel('speed (m/s)'); title(rd_names{s}); grid on;
    end
end

%% ── 图3a：SCR 时序（单独可视化）────────────────────────────────────────────
figure('Name','SCR Timing Sequence');
plot(t_frames,SCR_A,'k-','LineWidth',1.1,'DisplayName','A MTI Filtering'); hold on;
plot(t_frames,SCR_B,'g--','LineWidth',1.1,'DisplayName','B Wide Notch');
plot(t_frames,SCR_C,'b-','LineWidth',1.4,'DisplayName','C RD Adaptive Soft Notch');
plot(t_frames,SCR_D,'r-','LineWidth',1.7,'DisplayName','D Time-Frequency Ridge');
plot(t_frames,SCR_E,'m-','LineWidth',1.2,'DisplayName','E Subspace');
plot(t_frames,SCR_F,'c--','LineWidth',1.2,'DisplayName','F No-IMU Ridge');
xlabel('times (s)'); ylabel('SCR (dB)');
title(sprintf('Signal-to-Interference Ratio Timing of Each Frame（Mean: A=%.1f, B=%.1f, C=%.1f, D=%.1f dB）', ...
    mean(SCR_A), mean(SCR_B), mean(SCR_C), mean(SCR_D)));
legend('Location','best'); grid on;

%% ── 图3b：TERI 时序（单独可视化）──────────────────────────────────────────
figure('Name','TERI Timing Sequence','Position',[40,80,800,420]);
plot(t_frames,TERI_A,'k-','LineWidth',0.8,'DisplayName','A MTI Filtering'); hold on;
plot(t_frames,TERI_B,'g--','LineWidth',0.8,'DisplayName','B Wide Notch');
plot(t_frames,TERI_C,'b-','LineWidth',1.4,'DisplayName','C RD Adaptive Soft Notch');
plot(t_frames,TERI_D,'r-','LineWidth',1.7,'DisplayName','D Time-Frequency Ridge');
plot(t_frames,TERI_E,'m-','LineWidth',1.0,'DisplayName','E Subspace');
plot(t_frames,TERI_F,'c--','LineWidth',1.0,'DisplayName','F No-IMU Ridge');
yline(1.0,'k--','Benchmark','LineWidth',1.0);
xlabel('times (s)'); ylabel('Target Energy Retention Rate TERI');
title(sprintf('Target Energy Retention Rate（Mean: A=%.3f, B=%.3f, C=%.3f, D=%.3f）', ...
    mean(TERI_A), mean(TERI_B), mean(TERI_C), mean(TERI_D)));
legend('Location','best'); grid on;

%% ── 图3c：数据驱动杂波总展宽 σ_c 综合分析（单独可视化）────────────────────
figure('Name','Clutter Total Bandwidth Analysis');

plot(t_frames, sigma_c_total_arr,'m-','LineWidth',1.6,'DisplayName','σ_cData-driven(EMA)'); hold on;
plot(t_frames, clutter_model.sigma_c_est,'b--','LineWidth',1.2,'DisplayName','σ_c prior lower bound (platform)');
plot(t_frames, sigma_c_all,'g-','LineWidth',1.0,'DisplayName','σ_c final');
xlabel('times (s)'); ylabel('σ_c (m/s)');
title('σ_c Timing Comparison: Data-Driven EMA vs Prior Lower Bound vs Final Fusion ...');
legend('Location','best'); grid on;



% sgtitle('数据驱动杂波总展宽 σ_c 综合分析（时序+分布）','FontSize',12,'FontWeight','bold');

%% ── 图4：四方案性能柱状对比（消融）────────────────────────────────────────
methods = {'A-MTI Filtering','B Wide Notch','C-RD Adaptive Soft Notch','D-Time-Frequency Ridge','E Subspace (RPCA/eigencanceler)','F Data-driven Ridge (no IMU)'};
SCR_vals = [mean(SCR_A),mean(SCR_B),mean(SCR_C),mean(SCR_D),mean(SCR_E),mean(SCR_F)];
FAR_vals = [FAR_A,FAR_B,FAR_C,FAR_D,FAR_E,FAR_F]*1e4;
vmask_f4 = peak_found;
if sum(vmask_f4) > 0
    TERI_vals = [mean(TERI_A(vmask_f4)), mean(TERI_B(vmask_f4)), ...
        mean(TERI_C(vmask_f4)), mean(TERI_D(vmask_f4)), mean(TERI_E(vmask_f4)), mean(TERI_F(vmask_f4))];
else
    TERI_vals = [mean(TERI_A), mean(TERI_B), mean(TERI_C), mean(TERI_D), mean(TERI_E), mean(TERI_F)];
end
figure('Name','Average SCR (dB)');

bar(1:6,SCR_vals,0.5,'FaceColor',[0.3 0.6 0.9]);
set(gca,'XTick',1:6,'XTickLabel',methods,'XTickLabelRotation',12);
ylabel('Average SCR (dB)'); title('Confidence Correspondence Comparison'); grid on;
for i=1:6; text(i,SCR_vals(i),sprintf('%.1f',SCR_vals(i)),'Hor','center','Vert','bottom','FontSize',9); end

figure('Name','FAR (×10^{-4})');

bar(1:6,FAR_vals,0.5,'FaceColor',[0.9 0.4 0.3]);
set(gca,'XTick',1:6,'XTickLabel',methods,'XTickLabelRotation',12);
ylabel('FAR (×10^{-4})'); title('False Alarm Rate Comparison'); grid on;
for i=1:6; text(i,FAR_vals(i),sprintf('%.2f',FAR_vals(i)),'Hor','center','Vert','bottom','FontSize',9); end

figure('Name','TERI');

bar(1:6,TERI_vals,0.5,'FaceColor',[0.2 0.7 0.4]);
set(gca,'XTick',1:6,'XTickLabel',methods,'XTickLabelRotation',12);
ylim([0, max(1.2, max(TERI_vals)*1.15)]);
ylabel('TERI'); title('Target Energy Retention Rate'); grid on;
for i=1:6; text(i,TERI_vals(i),sprintf('%.3f',TERI_vals(i)),'Hor','center','Vert','bottom','FontSize',9); end
yline(1.0,'k--','LineWidth',1.0);
% sgtitle('四方案性能对比','FontSize',12,'FontWeight','bold');

%% ── 图5a~5f：方案 D 处理流水线各步骤分别展示 ──────────────────────────
if isfield(MID,'diagD') && isfield(MID.diagD,'Ptf') && isfield(MID.diagD,'G')
    dg = MID.diagD;
    tf_t   = ((1:dg.Nf)-1)*tfc.H*PRI*1e3;
    Ptf_dB = 10*log10(dg.Ptf/(max(dg.Ptf(:))+eps)+eps);
    Ptf_masked = dg.Ptf .* dg.G.^2;
    Pmsk_dB= 10*log10(Ptf_masked/(max(dg.Ptf(:))+eps)+eps);

    % ① 输入 RD（MOCO 后，未滤波）
    figure('Name','D-① 输入RD');
    imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_pre)+eps)); axis xy; colorbar; clim([-10,100]); hold on;
    if isfield(MID,'d_ref_k_arr')
        for kk = 1:MID.N_tgt
            plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)),'ro','MarkerSize',6,'LineWidth',1);
        end
    else
        plot(range_axis(MID.r_ref), v_tf(MID.d_ref),'r^','MarkerSize',10,'LineWidth',1.8);
    end
    xlabel('距离 (m)'); ylabel('速度 (m/s)'); title('① 输入 RD（MOCO 后，未滤波）'); grid on;

    % ② 慢时间 STFT 时频谱（杂波距离区）
    figure('Name','D-② STFT时频谱');
    imagesc(tf_t, v_tf, Ptf_dB); axis xy; colorbar; clim([-40,0]); colormap(gca,parula);
    xlabel('CPI 内时间 (ms)'); ylabel('速度 (m/s)');
    title('② 慢时间 STFT 时频谱（杂波距离区）'); grid on;

    % ③ IMU先验 + Viterbi 脊线追踪
    figure('Name','D-③ Viterbi脊线追踪');
    imagesc(tf_t, v_tf, Ptf_dB); axis xy; colorbar; clim([-40,0]); colormap(gca,parula); hold on;
    plot(tf_t, v_tf(dg.k_pred),'b--','LineWidth',1.5,'DisplayName','IMU先验');
    plot(tf_t, dg.v_ridge,'r-','LineWidth',2.0,'DisplayName','Viterbi脊线');
    plot(tf_t, dg.v_ridge+dg.B_flat_v,'m:','LineWidth',1.0,'DisplayName','平顶边界');
    plot(tf_t, dg.v_ridge-dg.B_flat_v,'m:','LineWidth',1.0);
    xlabel('CPI 内时间 (ms)'); ylabel('速度 (m/s)');
    title('③ IMU先验 + Viterbi 脊线追踪'); legend('Location','best'); grid on;

    % ④ 目标保持型时频掩模 G
    figure('Name','D-④ 时频掩模G');
    imagesc(tf_t, v_tf, dg.G); axis xy; colorbar; clim([0,1]); colormap(gca, flipud(gray)); hold on;
    plot(tf_t, dg.v_ridge,'r-','LineWidth',1.2);
    xlabel('CPI 内时间 (ms)'); ylabel('速度 (m/s)');
    title('④ 目标保持型时频掩模 G（白=保留,黑=抑制）'); grid on;

    % ⑤ 掩模后时频谱（杂波脊线被压制）
    figure('Name','D-⑤ 掩模后时频谱');
    imagesc(tf_t, v_tf, Pmsk_dB); axis xy; colorbar; clim([-40,0]); colormap(gca,parula);
    xlabel('CPI 内时间 (ms)'); ylabel('速度 (m/s)');
    title('⑤ 掩模后时频谱（杂波脊线被压制）'); grid on;

    % ⑥ 输出 RD（两次滤波后）
    figure('Name','D-⑥ 输出RD');
    imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_D)+eps)); axis xy; colorbar; clim([-10,100]); hold on;
    if isfield(MID,'d_ref_k_arr')
        for kk = 1:MID.N_tgt
            plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)),'ro','MarkerSize',6,'LineWidth',1,'DisplayName',sprintf('目标%d',kk));
        end
    else
        plot(range_axis(MID.r_ref), v_tf(MID.d_ref),'r^','MarkerSize',10,'LineWidth',1.8,'DisplayName','目标');
    end

    if isfield(MID,'protect_cells_D') && ~isempty(MID.protect_cells_D)
        plot(range_axis(MID.protect_cells_D(:,2)), v_tf(MID.protect_cells_D(:,1)), ...
            'gs','MarkerSize',7,'LineWidth',1.2,'DisplayName','速度门控后候选');
    end
    xlabel('距离 (m)'); ylabel('速度 (m/s)'); title('⑥ 输出 RD（两次滤波后）'); legend('Location','best'); grid on;
end

%% ── 图6a~6b：时频面目标标注分别展示 ──────────────────────────────────
if isfield(MID,'diagD') && isfield(MID.diagD,'Ptf_tgt')
    dg = MID.diagD;
    tf_t = ((1:dg.Nf)-1)*tfc.H*PRI*1e3;

    % (a) 杂波距离区时频谱 + 目标速度参考线
    figure('Name','杂波距离区时频谱');
    imagesc(tf_t, v_tf, 10*log10(dg.Ptf/(max(dg.Ptf(:))+eps)+eps)); axis xy; colorbar;
    clim([-40,0]); colormap(gca,parula); hold on;
    plot(tf_t, dg.v_ridge,'r-','LineWidth',1.8,'DisplayName','杂波脊线');
    yline(MID.v_tgt,'g--',sprintf('目标多普勒 %.1f m/s',MID.v_tgt),'LineWidth',1.6,'FontSize',9);
    xlabel('CPI 内时间 (ms)'); ylabel('速度 (m/s)');
    title('杂波距离区时频谱（目标不在此距离，仅标其多普勒）'); legend('Location','best'); grid on;

    % (b) 目标距离区时频谱
    figure('Name','目标距离区时频谱');
    Pt_dB = 10*log10(dg.Ptf_tgt/(max(dg.Ptf_tgt(:))+eps)+eps);
    imagesc(tf_t, v_tf, Pt_dB); axis xy; colorbar; clim([-40,0]); colormap(gca,parula); hold on;
    yline(MID.v_tgt,'r-',sprintf('目标 %.1f m/s',MID.v_tgt),'LineWidth',1.8,'FontSize',9, ...
        'LabelHorizontalAlignment','left');
    plot(tf_t, dg.v_ridge,'w:','LineWidth',1.2,'DisplayName','杂波脊线(参考)');
    xlabel('CPI 内时间 (ms)'); ylabel('速度 (m/s)');
    title(sprintf('目标距离区时频谱（R≈%.0fm）：目标为稳定水平脊线', MID.tgt_R));
    legend('Location','best'); grid on;
end

% %% ── 图7：距离自适应 CNR 权重 w(r) 可视化 ──────────────────────────────────
% %  （σ_c 展宽已由图3c独立展示，此处仅保留距离权重）
% if isfield(MID,'diagD') && isfield(MID.diagD,'w_r')
%     dg = MID.diagD;
%     figure('Name','距离自适应权重','Position',[40,80,760,360]);
%     plot(range_axis, dg.w_r,'b-','LineWidth',1.6); hold on;
%     xline(range_axis(min(end,MID.tgt_bin)),'r--',...
%         sprintf('目标 R=%.0fm',MID.tgt_R),'LineWidth',1.2,'FontSize',8);
%     ylim([0,1.05]); xlabel('距离 (m)'); ylabel('杂波存在权重 w(r)');
%     title(sprintf('距离自适应 CNR 权重 w(r)（第%d帧）：w≈0=无杂波不放陷波; w≈1=强杂波全压制', f_mid));
%     grid on;
% end

%% ── 图8：两次滤波流水线对比可视化（含目标真值标注）─────────────────────────
% % 此图直观呈现"检测发生在哪个阶段"的差异：
% %   ① 输入RD（未滤波）→ 目标可见但杂波强
% %   ② 保守预压制（B_core < B_flat）→ 杂波核心被压，目标在边缘区仍可见
% %   ③ 完整TF掩模后（恢复前）→ 杂波强力压制，但目标也可能被压
% %   ④ 最终输出（恢复后）→ 杂波已压，目标能量由速度门控候选驱动恢复
% %  所有面板均标注目标真值位置（红色三角），便于对比各阶段目标保持情况
% if isfield(MID,'rd_pass1') && isfield(MID,'rd_D_pure')
%     clim_rd = [-10, 100];
%     figure('Name','两次滤波流水线对比','Position',[40,80,1480,680]);
% 
%     subplot(2,2,1);
%     imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_pre)+eps)); axis xy; colorbar;
%     clim(clim_rd); hold on;
%     % 目标真值标注：红色三角标记所有目标的理论位置
%     if isfield(MID,'d_ref_k_arr') && isfield(MID,'r_ref_k_arr')
%         for kk = 1:MID.N_tgt
%             plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)), ...
%                 'r^','MarkerSize',10,'LineWidth',1.8, ...
%                 'DisplayName',sprintf('目标#%d真值',kk));
%         end
%     end
%     % 标注 B_flat 与 B_core 边界
%     yline(MID.mu_c + MID.B_flat, 'm--','LineWidth',1.4);
%     yline(MID.mu_c - MID.B_flat, 'm--','B_{flat}边界','LineWidth',1.4);
%     yline(MID.mu_c + MID.bw_3dB, 'c-','LineWidth',1.2);
%     yline(MID.mu_c - MID.bw_3dB, 'c-','B_{3dB}边界','LineWidth',1.2);
%     xlabel('距离 (m)'); ylabel('速度 (m/s)');
%     title(sprintf('① 输入RD（MOCO后，未滤波）\n品红虚线=B_{flat}=%.1fm/s  青线=B_{3dB}=%.1fm/s', ...
%         MID.B_flat, MID.bw_3dB));
%     legend('Location','best'); grid on;
% 
%     subplot(2,2,2);
%     imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_pass1)+eps)); axis xy; colorbar;
%     clim(clim_rd); hold on;
%     % 目标真值标注（同上）
%     if isfield(MID,'d_ref_k_arr') && isfield(MID,'r_ref_k_arr')
%         for kk = 1:MID.N_tgt
%             plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)), ...
%                 'r^','MarkerSize',10,'LineWidth',1.8);
%         end
%     end
%     % 速度门控边界
%     yline(MID.mu_c + (MID.B_flat_gate + DELTA_MARGIN), 'm--','LineWidth',1.4);
%     yline(MID.mu_c - (MID.B_flat_gate + DELTA_MARGIN), 'm--', 'LineWidth',1.4);
%     yline(MID.mu_c + MID.B_flat_used, 'y:','LineWidth',1.2);
%     yline(MID.mu_c - MID.B_flat_used, 'y:', 'LineWidth',1.2);
%     xlabel('距离 (m)'); ylabel('速度 (m/s)');
%     title(sprintf('② 第一次保守压制（3dB带宽=%.1fm/s）\n门控带已与TF掩模宽度对齐 → 带内杂波不再误判为目标', MID.bw_3dB));
%     legend('Location','best'); grid on;
% 
%     subplot(2,2,3);
%     imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_D_pure)+eps)); axis xy; colorbar;
%     clim(clim_rd); hold on;
%     % 目标真值标注：观察目标是否被TF掩模误压制
%     if isfield(MID,'d_ref_k_arr') && isfield(MID,'r_ref_k_arr')
%         for kk = 1:MID.N_tgt
%             plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)), ...
%                 'r^','MarkerSize',10,'LineWidth',1.8, ...
%                 'DisplayName',sprintf('目标#%d',kk));
%         end
%     end
%     xlabel('距离 (m)'); ylabel('速度 (m/s)');
%     title(sprintf('③ 完整TF掩模后（恢复前）\n目标可能被压制（B_{flat}=%.1fm/s）', MID.B_flat));
%     legend('Location','best'); grid on;
% 
%     subplot(2,2,4);
%     imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_D)+eps)); axis xy; colorbar;
%     clim(clim_rd); hold on;
%     % 目标真值标注：验证最终输出后目标是否被成功恢复
%     if isfield(MID,'d_ref_k_arr') && isfield(MID,'r_ref_k_arr')
%         for kk = 1:MID.N_tgt
%             plot(range_axis(MID.r_ref_k_arr(kk)), v_tf(MID.d_ref_k_arr(kk)), ...
%                 'r^','MarkerSize',10,'LineWidth',1.8, ...
%                 'DisplayName',sprintf('目标#%d真值',kk));
%         end
%     end
%     if isfield(MID,'protect_cells_D') && ~isempty(MID.protect_cells_D)
%         plot(range_axis(MID.protect_cells_D(:,2)), v_tf(MID.protect_cells_D(:,1)), ...
%             'gs','MarkerSize',8,'LineWidth',1.5,'DisplayName','速度门控后候选（已恢复能量）');
%     end
%     xlabel('距离 (m)'); ylabel('速度 (m/s)');
%     title(sprintf('④ 最终输出（精准恢复后）\n候选: 原始%d → 速度门控后%d', ...
%         MID.n_cand_raw, MID.n_cand_kept));
%     legend('Location','best'); grid on;
% 
%     sgtitle(sprintf(['两次滤波流水线（第%d帧）—— 红色三角=目标真值位置\n'...
%         '核心：在保守压制结果上检测候选 + 速度域门控去虚警 + 全陷波精准恢复'], f_mid), ...
%         'FontSize',11,'FontWeight','bold');
% end

%% ── 图9：检测性能三合一（Pd-Pfa / FAR-Pfa / ROC）────────────────────────
% figure('Name','检测性能分析','Position',[40,80,1000,520]);
method_names = {'A-MTI Filtering','B Wide Notch','C-RD Soft Notch Filter','D-Time-Frequency Ridge','E Subspace (RPCA/eigencanceler)','F Data-driven Ridge (no IMU)'};
line_sty = {'k-o','g-s','b-^','r-d','m-p','c-x'};
% line_sty = {'k-','g-','b-','r-'};
x_pfa = PFA_SWEEP;

% ── 子图1：Pd vs Pfa（半对数）──
figure('Name',' Pd vs Pfa');
for s_idx = 1:6
    semilogx(x_pfa, Pd_roc(:,s_idx), line_sty{s_idx}, ...
        'LineWidth',1.6,'MarkerSize',5,'MarkerFaceColor','auto', ...
        'DisplayName',method_names{s_idx}); hold on;
end
xline(Pfa_work, 'k:', 'Pfa_{nom}=10^{-4}', 'HandleVisibility','off');
xlabel('Pfa'); ylabel('Detection Probability P_d'); title('(a) Pd vs Pfa');
legend('Location','southeast'); grid on; xlim([min(x_pfa)/2, max(x_pfa)*2]); ylim([0,1.05]);

% % ── 子图2：FAR vs Pfa（双对数）──
% subplot(1,3,2);
figure('Name','  FAR Curve');
for s_idx = 1:6
    loglog(x_pfa, FAR_roc(:,s_idx), line_sty{s_idx}, ...
        'LineWidth',1.6,'MarkerSize',5,'MarkerFaceColor','auto'); hold on;
end
loglog(x_pfa, x_pfa, 'k:', 'LineWidth',1.0, 'HandleVisibility','off');  % 对角线=理想映射
xlabel('Pfa'); ylabel('FAR'); title('(b) FAR vs Pfa');
grid on; xlim([min(x_pfa)/2, max(x_pfa)*2]); ylim([min(FAR_roc(FAR_roc>0))/2, 1]);
legend('Location','southeast'); grid on; xlim([min(x_pfa)/2, max(x_pfa)*2]); ylim([0,1.05]);
% ── 子图3：ROC 曲线（Pd vs FAR）──
figure('Name','ROC Curve');
for s_idx = 1:6
    loglog(FAR_roc(:,s_idx), Pd_roc(:,s_idx), line_sty{s_idx}, ...
        'LineWidth',1.8,'MarkerSize',5,'MarkerFaceColor','auto', ...
        'DisplayName',method_names{s_idx}); hold on;
end
for s_idx = 1:6
    loglog(FAR_roc(pfa_ref_idx,s_idx), Pd_roc(pfa_ref_idx,s_idx), ...
        'ks','MarkerSize',11,'LineWidth',2.5);
end
% FAR=10^{-4} 参考线和 Pd=0.9 参考线
xline(1e-4, 'k:', 'FAR=10^{-4}', 'HandleVisibility','off');
yline(0.9,  'k:', 'P_d=0.9',    'HandleVisibility','off');
xlabel(' FAR'); ylabel('Detection Probability P_d'); title('(c) ROC Cruve');
legend(method_names, 'Location','southeast');
grid on;
xlim([min(FAR_roc(FAR_roc>0))/2, 1]);
ylim([0, 1.05]);

% sgtitle(sprintf('检测性能综合分析（AUC: A=%.3f B=%.3f C=%.3f D=%.3f）', AUC), ...
%     'FontSize',12,'FontWeight','bold');

%% ── 图12a~12c：多目标检测性能分别展示 ────────────────────────────────
tgt_speeds = arrayfun(@(k) norm(tgt(k).vt), 1:N_tgt);
if nv > 5 && exist('SCR_D_tgt','var')
    % (a) 逐目标 SCR_D
    figure('Name','逐目标 SCR_D');
    bar(1:N_tgt, mean(SCR_D_tgt(vmask,:),1), 0.6, 'FaceColor',[0.9 0.3 0.2]); hold on;
    errorbar(1:N_tgt, mean(SCR_D_tgt(vmask,:),1), std(SCR_D_tgt(vmask,:),[],1), 'k.','LineWidth',1.2);
    set(gca,'XTick',1:N_tgt); xlabel('Target #'); ylabel('Mean SCR_D (dB)');
    title('Per-target SCR_D (Scheme D, valid frames)'); grid on;
    xlim([0.5, N_tgt+0.5]);

    % (b) 逐目标 TERI_D
    figure('Name','逐目标 TERI_D');
    bar(1:N_tgt, mean(TERI_D_tgt(vmask,:),1), 0.6, 'FaceColor',[0.2 0.7 0.4]); hold on;
    yline(1.0,'k--','LineWidth',1.0);
    set(gca,'XTick',1:N_tgt); xlabel('Target #'); ylabel('Mean TERI_D');
    title('Per-target TERI_D (Scheme D, valid frames)'); grid on;
    xlim([0.5, N_tgt+0.5]);

    % (c) SCR_D vs 目标速度
    figure('Name','SCR_D vs 目标速度');
    scatter(tgt_speeds, mean(SCR_D_tgt(vmask,:),1), 80, 'filled');
    xlabel('Target Speed (m/s)'); ylabel('Mean SCR_D (dB)');
    title('SCR_D vs Target Speed'); grid on;
    for kk = 1:N_tgt
        text(tgt_speeds(kk)+0.1, mean(SCR_D_tgt(vmask,kk),1), sprintf('#%d',kk),'FontSize',9);
    end
end

% %% ── 图13：慢速目标多普勒剖面对比（核心验证：D vs C 的保持优势）─────────────
% if isfield(MID,'rd_A') && isfield(MID,'rd_D') && isfield(MID,'rd_pre')
%     mu_c_mid = MID.mu_c;
%     v_rel = abs(abs(MID.tgt_vr_arr) - abs(mu_c_mid));
%     v_rel(v_rel < 0.3) = inf;
%     [~, slowest_idx] = min(v_rel);
%     r_tgt = MID.r_ref_k_arr(slowest_idx);
%     figure('Name','慢速目标多普勒剖面对比','Position',[40,80,1000,520]);
%     dop_axis = v_tf(:);
%     p_max = max(abs(MID.rd_pre(:, r_tgt)).^2) + eps;
%     p_pre_dB = 10*log10(abs(MID.rd_pre(:, r_tgt)).^2 / p_max + eps);
%     p_A_dB   = 10*log10(abs(MID.rd_A(:, r_tgt)).^2 / p_max + eps);
%     p_B_dB   = 10*log10(abs(MID.rd_B(:, r_tgt)).^2 / p_max + eps);
%     p_C_dB   = 10*log10(abs(MID.rd_C(:, r_tgt)).^2 / p_max + eps);
%     p_D_dB   = 10*log10(abs(MID.rd_D(:, r_tgt)).^2 / p_max + eps);
%     subplot(1,2,1);
%     plot(dop_axis, p_pre_dB, 'Color',[0.4 0.4 0.4],'LineWidth',1.0,'DisplayName','原始RD'); hold on;
%     plot(dop_axis, p_A_dB, 'c-','LineWidth',1.0,'DisplayName','A MTI');
%     plot(dop_axis, p_B_dB, 'Color',[0 0.6 0],'LineWidth',1.0,'DisplayName','B宽陷波');
%     plot(dop_axis, p_C_dB, 'b-','LineWidth',1.6,'DisplayName','C RD软陷波');
%     plot(dop_axis, p_D_dB, 'r-','LineWidth',2.0,'DisplayName','D时频脊线');
%     v_tgt_mid = MID.tgt_vr_arr(slowest_idx);
%     xline(v_tgt_mid,'r--',sprintf('T#%d v_{rad}=%.1f',slowest_idx,v_tgt_mid),'LineWidth',1.2,'FontSize',8,'LabelHorizontalAlignment','left');
%     xline(mu_c_mid,'m:','\mu_c','LineWidth',1.2,'FontSize',8);
%     xline(mu_c_mid+MID.B_flat,'m--','+B_{flat}','LineWidth',0.8,'FontSize',7);
%     xline(mu_c_mid-MID.B_flat,'m--','-B_{flat}','LineWidth',0.8,'FontSize',7);
%     xlabel('速度(m/s)'); ylabel('归一化功率(dB)');
%     title(sprintf('多普勒@R=%.0fm(帧%d)',range_axis(r_tgt),f_mid));
%     legend('Location','southwest','FontSize',7); grid on; xlim([min(v_tf),max(v_tf)]); ylim([-60,5]);
%     subplot(1,2,2);
%     zoom_margin = max(abs(v_tgt_mid-mu_c_mid)*2.5,3.0);
%     zm = dop_axis >= mu_c_mid-zoom_margin & dop_axis <= mu_c_mid+zoom_margin;
%     plot(dop_axis(zm),p_pre_dB(zm),'Color',[0.4 0.4 0.4],'LineWidth',1.0); hold on;
%     plot(dop_axis(zm),p_C_dB(zm),'b-','LineWidth',1.8);
%     plot(dop_axis(zm),p_D_dB(zm),'r-','LineWidth',2.5);
%     xline(v_tgt_mid,'r--',sprintf('T#%d',slowest_idx),'LineWidth',1.2,'FontSize',8,'LabelHorizontalAlignment','left');
%     xline(mu_c_mid,'m:','\mu_c','LineWidth',1.2,'FontSize',8);
%     [~,tgt_bin_v] = min(abs(dop_axis-v_tgt_mid));
%     delta_dB = p_D_dB(tgt_bin_v)-p_C_dB(tgt_bin_v);
%     text(v_tgt_mid+0.3,p_D_dB(tgt_bin_v)+1,sprintf('D-C=%.1fdB',delta_dB),'Color','r','FontWeight','bold');
%     xlabel('速度(m/s)'); ylabel('dB');
%     title(sprintf('局部放大 D-C=%.1fdB@目标',delta_dB));
%     grid on; xlim(mu_c_mid+[-1,1]*zoom_margin); ylim([-50,5]);
%     legend({'原始RD','C RD软陷波','D时频脊线'},'Location','southwest','FontSize',7);
%     sgtitle(sprintf('慢速目标多普勒剖面：D保持目标能量(红线) vs C削弱目标(蓝线)'),'FontSize',11,'FontWeight','bold');
% end

%% ── 图15：慢速目标检测概率(Pd_slow)与快速目标检测概率(Pd_fast)对比 ────────
%  展示四方案在杂波带内目标(Pd_slow)和杂波带外目标(Pd_fast)上的检测性能差异，
%  重点突出方案D对慢速目标的保持优势。
method_names_4 = {'A-MTI Filtering','B Wide Notch','C-RD Soft Notch Filter','D-Time-Frequency Ridge','E Subspace','F No-IMU Ridge'};
figure('Name','Slow Target Detection Probability Pd_slow vs Pd_fast');
bar_width = 0.35;
x_pos = 1:6;
% 并排柱状图：Pd_slow 和 Pd_fast 在每个方案下并列
bar(x_pos - bar_width/2, Pd_slow, bar_width, 'FaceColor',[0.8 0.3 0.3]); hold on;
bar(x_pos + bar_width/2, Pd_fast, bar_width, 'FaceColor',[0.3 0.6 0.8]);
% 标注数值
for i = 1:6
    text(x_pos(i) - bar_width/2, Pd_slow(i) + 0.015, sprintf('%.3f',Pd_slow(i)), ...
        'Hor','center','Vert','bottom','FontSize',9,'Color',[0.8 0.3 0.3]);
    text(x_pos(i) + bar_width/2, Pd_fast(i) + 0.015, sprintf('%.3f',Pd_fast(i)), ...
        'Hor','center','Vert','bottom','FontSize',9,'Color',[0.3 0.6 0.8]);
end
% D vs C 改善标注
dc_improve = (Pd_slow(4) - Pd_slow(3)) * 100;
plot([3 4], [Pd_slow(3) Pd_slow(4)]*1.05, 'k-', 'LineWidth',1.2);
text(3.5, max(Pd_slow(3:4))*1.08, sprintf('Method D has improved\n+%.1f%% (相对+%.0f%%)', ...
    dc_improve, (Pd_slow(4)/max(Pd_slow(3),eps)-1)*100), ...
    'Hor','center','FontSize',9,'Color','r','FontWeight','bold');
set(gca,'XTick',1:6,'XTickLabel',method_names_4,'XTickLabelRotation',12);
ylabel('Detection Probability P_d'); title('Comparison of Detection Probabilities between Slow and Fast Target Detection');
legend({'Pd_{slow} ', 'Pd_{fast} '}, 'Location','northeast');
grid on; ylim([0, max([Pd_slow,Pd_fast])*1.3+0.05]);

%% ── 图14：目标保持率 vs 目标-杂波间距（定量证明D对慢速目标的保持优势）─────
if exist('TERI_D_tgt','var') && exist('TERI_C_tgt','var') && exist('tgt_vr_all','var') && exist('mu_c_all','var')
    x_all=[]; y_all=[]; tgt_id_all=[];
    for f=1:num_frames
        for k=1:N_tgt
            dv = abs(tgt_vr_all(f,k)-mu_c_all(f));
            ratio = TERI_D_tgt(f,k)/max(TERI_C_tgt(f,k),eps);
            if isfinite(ratio) && ratio>0 && ratio<10 && dv<10
                x_all=[x_all;dv]; y_all=[y_all;ratio]; tgt_id_all=[tgt_id_all;k];
            end
        end
    end
    if length(x_all)>10
%         figure('Name','目标保持率与杂波间距关系');
% %         subplot(1,2,1); 
%         colors_tgt=lines(N_tgt); hold on;
%         for k=1:N_tgt
%             mk=tgt_id_all==k;
%             if sum(mk)>0; scatter(x_all(mk),y_all(mk),12,colors_tgt(k,:),'filled','DisplayName',sprintf('T%d(spd=%.1f)',k,tgt_speeds(k))); end
%         end
%         yline(1.0,'k--','TERI_C=TERI_D','LineWidth',1.2);
%         [xs,si]=sort(x_all); ys=y_all(si);
%         plot(xs,movmean(ys,max(3,round(length(xs)/20))),'r-','LineWidth',2.0,'DisplayName','趋势');
%         xlabel('|v_{rad} - \mu_c| (m/s)'); ylabel('TERI_D / TERI_C');
%         title('D/C保持率 vs 目标-杂波间距'); legend('Location','best','FontSize',7); grid on;
        figure()
        mr=zeros(N_tgt,1); 
        sr=mr;
        for k=1:N_tgt; mk=tgt_id_all==k;
            if sum(mk)>0; mr(k)=mean(y_all(mk)); sr(k)=std(y_all(mk)); end
        end
        bar(1:N_tgt,mr,0.6,'FaceColor',[0.8 0.3 0.3]); hold on;
        errorbar(1:N_tgt,mr,sr,'k.','LineWidth',1.2);
        yline(1.0,'k--','Benchmark','LineWidth',1.2);
        ylabel('TERI_D/TERI_C'); 
        xlabel('Target'); 
        title('Target-by-target D/C retention ratio');
        set(gca,'XTick',1:N_tgt,'XTickLabel',arrayfun(@(k)sprintf('T%d\n%.1f',k,tgt_speeds(k)),1:N_tgt,'UniformOutput',false));
        grid on; ylim([0,max(2,max(mr)*1.2)]);
        title('For slow targets , D is significantly better than C; for fast targets , the two are converging.','FontSize',12,'FontWeight','bold');
    end
end

%% =========================================================================
%  PART 12：保存实验数值结果与可视化图表（分目录保存）
% =========================================================================
%  设计说明：
%    - data/  : 数值实验结果（.mat 文件 + metrics_summary.txt 文本汇总）
%    - figs/  : 可视化图表（.png 用于浏览，.fig 用于 MATLAB 再编辑）
%    - 根目录 : index.html 图库索引
% =========================================================================
results_dir = sprintf('results_%s', datestr(now, 'yyyymmdd_HHMMSS'));
mkdir(results_dir);
mkdir(fullfile(results_dir, 'data'));
mkdir(fullfile(results_dir, 'figs'));

fprintf('\n[保存] 保存实验数值结果至 %s/data/ ...\n', results_dir);

% ── 12a. 保存数值实验结果（.mat 文件，分类存储）─────────────────────────

%  信杂比与能量保持率
save(fullfile(results_dir, 'data', 'metrics_SCR.mat'), ...
    'SCR_A','SCR_B','SCR_C','SCR_D','SCR_D_tgt', ...
    'TERI_A','TERI_B','TERI_C','TERI_D','TERI_D_tgt','TERI_C_tgt');

%  检测性能（ROC/AUC/Pd/FAR）
save(fullfile(results_dir, 'data', 'metrics_detection.mat'), ...
    'Pd_roc','FAR_roc','PFA_SWEEP','AUC', ...
    'Pd','FAR_A','FAR_B','FAR_C','FAR_D', ...
    'Pd_sep_roc','Pd_clut_roc','n_sep_frames','n_clut_frames', ...
    'Pd_slow','Pd_fast','pd_slow_cnt','pd_fast_cnt', ...
    'SCR_D0','TERI_D0','Pd_slow_D0','Pd_fast_D0', ...
    'SCR_Dfx','TERI_Dfx','Pd_slow_Dfx','Pd_fast_Dfx', ...
    'SCR_Dfn','TERI_Dfn','Pd_slow_Dfn','Pd_fast_Dfn', ...
    'FAR_D0','FAR_Dfx','FAR_Dfn','SCR_Nc','Pd_slow_Nc','Nc_list', ...
    'n_slow_tgtframes','n_fast_tgtframes');

%  杂波参数（展宽/中心/3dB带宽）
save(fullfile(results_dir, 'data', 'metrics_clutter.mat'), ...
    'sigma_c_total_arr','sigma_c_all','mu_c_all', ...
    'bw_3dB_arr','bw_3dB_ratio_arr','sigma_c_ema', ...
    'clutter_model','BW3_C','BW3_D');

%  目标与系统参数
save(fullfile(results_dir, 'data', 'params_target.mat'), ...
    'N_tgt','tgt','tgt_speeds','tgt_vr_all');
save(fullfile(results_dir, 'data', 'params_system.mat'), ...
    'Ob_time','fc','B','PRI','lambda','N_slow','N_fast','Num_chirp', ...
    'range_axis','v_tf','range_res','num_frames','step');

%  中间帧诊断数据（若有）
if exist('MID','var') && ~isempty(fieldnames(MID))
    save(fullfile(results_dir, 'data', 'mid_frame_diagnostic.mat'), 'MID');
end

% ── 12b. 生成文本格式指标汇总 ──────────────────────────────────────────
fid_summary = fopen(fullfile(results_dir, 'data', 'metrics_summary.txt'), 'w');
fprintf(fid_summary, '=============================================================\n');
fprintf(fid_summary, '  时频脊线追踪 + 目标保持型时频掩模 — 仿真实验指标汇总\n');
fprintf(fid_summary, '  生成时间: %s\n', datestr(now));
fprintf(fid_summary, '=============================================================\n\n');

fprintf(fid_summary, '【系统参数】\n');
fprintf(fid_summary, '  载频 fc=%.1f GHz | 带宽 B=%.1f MHz | PRI=%.1f us\n', ...
    fc/1e9, B/1e6, PRI*1e6);
fprintf(fid_summary, '  波长 lambda=%.4f m | 距离分辨率=%.2f m\n', lambda, range_res);
fprintf(fid_summary, '  CPI 脉冲数 N_slow=%d | 总脉冲数=%d | 帧数=%d\n', ...
    N_slow, Num_chirp, num_frames);
fprintf(fid_summary, '  目标数=%d | 速度范围 [%.2f, %.2f] m/s\n\n', ...
    N_tgt, min(tgt_speeds), max(tgt_speeds));

fprintf(fid_summary, '【四方案性能对比】\n');
fprintf(fid_summary, '%-26s %9s %9s %9s %9s %9s %9s\n', ...
    '指标','A-MTI','B宽陷波','C-RD软陷波','D-时频脊线','E-子空间','F-无IMU脊线');
fprintf(fid_summary, '%s\n', repmat('-',1,84));
fprintf(fid_summary, '%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','平均SCR (dB)', ...
    mean(SCR_A), mean(SCR_B), mean(SCR_C), mean(SCR_D), mean(SCR_E), mean(SCR_F));
fprintf(fid_summary, '%-26s %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n','FAR', ...
    FAR_A,FAR_B,FAR_C,FAR_D,FAR_E,FAR_F);
fprintf(fid_summary, '%-26s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n','Pd', ...
    Pd(1),Pd(2),Pd(3),Pd(4),Pd(5),Pd(6));
fprintf(fid_summary, '%-26s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n','AUC (标准ROC)', AUC);
fprintf(fid_summary, '%-26s %9s %9.2f %9.2f %9.2f %9.2f %9.2f\n','SCR增益 vs A (dB)','-', ...
    mean(SCR_B)-mean(SCR_A), mean(SCR_C)-mean(SCR_A), mean(SCR_D)-mean(SCR_A), mean(SCR_E)-mean(SCR_A), mean(SCR_F)-mean(SCR_A));
fprintf(fid_summary, '%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','TERI', ...
    mean(TERI_A), mean(TERI_B), mean(TERI_C), mean(TERI_D), mean(TERI_E), mean(TERI_F));
fprintf(fid_summary, '%-26s %9s %9s %9.2f %9.2f %9s %9s\n','杂波-3dB谱宽(m/s)','-','-', ...
    mean(BW3_C(BW3_C>0)), mean(BW3_D(BW3_D>0)), '-', '-');
fprintf(fid_summary, '%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','Pd(slow,in-clut)', ...
    Pd_slow);
fprintf(fid_summary, '%-26s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n','Pd(fast,out-clut)', ...
    Pd_fast);
fprintf(fid_summary, '%s\n\n', repmat('=',1,66));

fprintf(fid_summary, '【关键对比 D vs C】\n');
fprintf(fid_summary, '  SCR 改善:        %+.2f dB\n', mean(SCR_D)-mean(SCR_C));
fprintf(fid_summary, '  TERI 比值(D/C):  %.3f\n', ...
    mean(TERI_D)/max(mean(TERI_C),eps));
fprintf(fid_summary, '  Slow Pd 改善:    %+.1f%%\n', (Pd_slow(4)-Pd_slow(3))*100);
fprintf(fid_summary, '  FAR 差异(D-C):   %+.2e\n', FAR_D-FAR_C);
fprintf(fid_summary, '  Pd 差异(D-C):    %+.2f\n', Pd(4)-Pd(3));
fclose(fid_summary);
fprintf('[保存] 指标汇总: data/metrics_summary.txt\n');

% ── 12c. 保存可视化图表（分图保存，标注命名）─────────────────────────────
fig_handles = findobj('Type', 'figure');
n_figs = length(fig_handles);
fprintf('\n[保存] 保存 %d 张可视化图表至 %s/figs/ ...\n', n_figs, results_dir);

% 图号 → 英文标注文件名映射表（避免中文字符在跨平台时的问题，
% 同时提供清晰的语义化命名用于论文插图引用）
fig_name_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
fig_name_map(1)  = 'fig1a_tf_spectrogram';             % 时频谱
fig_name_map(2)  = 'fig1b_tf_mask';                    % 时频掩模G
fig_name_map(3)  = 'fig2a_rd_orig';                    % 原始RD（MOCO后）
fig_name_map(4)  = 'fig2b_rd_mti';                     % A-MTI滤波
fig_name_map(5)  = 'fig2c_rd_widenotch';               % B-宽硬陷波
fig_name_map(6)  = 'fig2d_rd_softnotch';               % C-RD软陷波
fig_name_map(7)  = 'fig2e_rd_tfridge';                 % D-时频脊线
fig_name_map(8)  = 'fig3_scr_timeseries';              % SCR时序
fig_name_map(9)  = 'fig4_teri_timeseries';             % TERI时序
fig_name_map(10) = 'fig5_clutter_spread_analysis';     % sigma_c展宽分析
fig_name_map(11) = 'fig6a_ablation_scr';               % 消融-SCR
fig_name_map(12) = 'fig6b_ablation_far';               % 消融-FAR
fig_name_map(13) = 'fig6c_ablation_teri';              % 消融-TERI
fig_name_map(14) = 'fig7a_pipeline_input_rd';          % D-①输入RD
fig_name_map(15) = 'fig7b_pipeline_stft';              % D-②STFT时频谱
fig_name_map(16) = 'fig7c_pipeline_viterbi';           % D-③Viterbi脊线
fig_name_map(17) = 'fig7d_pipeline_mask_G';            % D-④时频掩模G
fig_name_map(18) = 'fig7e_pipeline_masked_spec';       % D-⑤掩模后时频谱
fig_name_map(19) = 'fig7f_pipeline_output_rd';         % D-⑥输出RD
fig_name_map(20) = 'fig8a_tf_clutter_zone';            % 杂波距离区时频谱
fig_name_map(21) = 'fig8b_tf_target_zone';             % 目标距离区时频谱
fig_name_map(22) = 'fig9a_pd_vs_pfa';                  % Pd vs 名义Pfa
fig_name_map(23) = 'fig9b_roc_curve';                  % ROC曲线
fig_name_map(24) = 'fig10a_per_target_scr';            % 逐目标SCR_D
fig_name_map(25) = 'fig10b_per_target_teri';           % 逐目标TERI_D
fig_name_map(26) = 'fig10c_scr_vs_speed';              % SCR_D vs 速度
fig_name_map(27) = 'fig11_teri_ratio_vs_dist';         % D/C保持率 vs 间距
fig_name_map(28) = 'fig12_per_target_dc_ratio';        % 逐目标D/C保持率
fig_name_map(29) = 'fig13_slow_fast_detection_pd';     % Pd_slow vs Pd_fast

fig_names_saved  = cell(n_figs, 1);
fig_annotations  = cell(n_figs, 1);

for h_idx = 1:n_figs
    h = fig_handles(h_idx);
    orig_name = get(h, 'Name');
    fig_num = h.Number;

    % 优先使用映射表名称（英文标注），否则使用清理后的 Name
    if isKey(fig_name_map, fig_num)
        fname = fig_name_map(fig_num);
    elseif ~isempty(orig_name)
        fname = regexprep(orig_name, '[^\w一-鿿_-]', '_');
        fname = regexprep(fname, '_+', '_');
        fname = regexprep(fname, '^_|_$', '');
    else
        fname = sprintf('figure_%d', fig_num);
    end

    fig_names_saved{h_idx} = fname;
    fig_annotations{h_idx} = orig_name;

    % 分别保存为 PNG（可视化浏览）和 FIG（MATLAB 可编辑）
    saveas(h, fullfile(results_dir, 'figs', [fname '.png']));
    saveas(h, fullfile(results_dir, 'figs', [fname '.fig']));
    fprintf('  [%d/%d] %s  (%s)\n', h_idx, n_figs, [fname '.png'], orig_name);
end

% ── 12d. 生成 HTML 图库（方便快速浏览全部图表）─────────────────────────
fid = fopen(fullfile(results_dir, 'index.html'), 'w');
fprintf(fid, '<!DOCTYPE html><html><head><title>Experiment Results</title>\n');
fprintf(fid, '<style>body{font-family:Arial,sans-serif;margin:20px;background:#f5f5f5}');
fprintf(fid, 'h1{color:#333}.gallery{display:flex;flex-wrap:wrap;gap:20px}');
fprintf(fid, '.card{background:white;border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,0.1);padding:15px;width:100%%;max-width:900px}');
fprintf(fid, '.card h2{font-size:16px;margin:0 0 10px 0;color:#555}');
fprintf(fid, '.card img{width:100%%;height:auto;border:1px solid #ddd;border-radius:4px}');
fprintf(fid, '</style></head><body>\n');
fprintf(fid, '<h1>多目标蜂群实验 — 时频脊线追踪 + 目标保持型掩模</h1>\n');
fprintf(fid, '<p>生成时间: %s | 目标数: %d | 速度范围: [%.2f, %.2f] m/s</p>\n', ...
    datestr(now), N_tgt, min(tgt_speeds), max(tgt_speeds));
fprintf(fid, '<p><a href="data/metrics_summary.txt">[浏览指标汇总]</a> | ');
fprintf(fid, '<a href="data/">[浏览数据文件]</a></p>\n');
fprintf(fid, '<div class="gallery">\n');
for h_idx = 1:n_figs
    fprintf(fid, '<div class="card">\n');
    fprintf(fid, '  <h2>%d. %s</h2>\n', h_idx, fig_annotations{h_idx});
    fprintf(fid, '  <p>文件名: %s.png</p>\n', fig_names_saved{h_idx});
    fprintf(fid, '  <a href="figs/%s.png"><img src="figs/%s.png" alt="%s"></a>\n', ...
        fig_names_saved{h_idx}, fig_names_saved{h_idx}, fig_annotations{h_idx});
    fprintf(fid, '  <p><a href="figs/%s.fig">[下载 FIG]</a></p>\n', fig_names_saved{h_idx});
    fprintf(fid, '</div>\n');
end
fprintf(fid, '</div></body></html>');
fclose(fid);

fprintf('\n[保存] HTML 图库: %s/index.html\n', results_dir);
fprintf('[保存] 数值数据: %s/data/ 目录下 (*.mat + summary.txt)\n', results_dir);
fprintf('[保存] 可视化图表: %s/figs/ 目录下 (*.png + *.fig)\n', results_dir);
fprintf('[保存] 实验数据与可视化图表已分目录保存完毕。\n\n');

fprintf('[完成] 所有实验数据与图表均已生成。\n');

%% =========================================================================
%  ====================  内 置 函 数 区  ===================================
%% =========================================================================

% ─────────────────────────────────────────────────────────────────────────
function rd_out = apply_conservative_notch(rd_in, v_axis, mu_c, B_core)
%APPLY_CONSERVATIVE_NOTCH  保守核心陷波（两次滤波第一阶段）
%  功能：对 RD 图施加一个窄陷波，仅压制杂波核心速度区 [μ_c-B_core, μ_c+B_core]，
%        保留边缘区（含 B_flat 附近的目标），供后续 CFAR 候选检测使用。
%  设计说明：
%    - 核心区增益 H_min=0.01（与主掩模 tf_ridge_filter 一致，确保核心杂波压制）
%    - 边缘区增益 = 1.0（完全不压制，目标在此区域保持原始能量）
%    - 陷波为矩形（硬边界），不做平滑过渡（第一阶段不需要精确边界）
H_min_p1 = 0.01;           % 核心区增益（-40 dB），与主掩模保持一致
v_col    = v_axis(:);       % 确保为列向量 [N_slow × 1]
H_vec    = ones(numel(v_col), 1);
core_idx = abs(v_col - mu_c) <= B_core;
H_vec(core_idx) = H_min_p1;
H        = repmat(H_vec, 1, size(rd_in, 2));   % [N_slow × N_fast]
rd_out   = rd_in .* H;
end

% ------------------------------------------------------------------------
function bw_3dB = compute_clutter_3dB_bw(rd, v_axis, ~, clut_bins)
%%COMPUTE_CLUTTER_3DB_BW  From RD domain, estimate the 3dB bandwidth of clutter peak
%  Input:
%    rd       - [N_slow x N_fast] RD data (after MOCO)
%    v_axis   - [N_slow x 1] velocity axis (m/s)
%    mu_c     - clutter center estimate (m/s)
%    clut_bins - range bins dominated by clutter
%  Output:
%    bw_3dB   - clutter 3dB bandwidth (m/s)
%
%  Principle: average Doppler spectrum over clut_bins -> normalize ->
%  search from peak to both sides for -3dB boundary
%  Replaces the fixed ratio BETA_CONSERVATIVE=0.8 with data-driven approach
v_axis = v_axis(:);
N_dop = numel(v_axis);

% Average Doppler spectrum over clutter range bins
cb = clut_bins(clut_bins>=1 & clut_bins<=size(rd,2));
if isempty(cb)
    cb = 1:size(rd,2);
end
pwr = mean(abs(rd(:, cb)).^2, 2);

% Remove noise floor and normalize
pwr = pwr - min(pwr);
pwr_norm = pwr / max(pwr);

% Search for continuous -3dB region around peak
[~, pk_idx] = max(pwr_norm);
lo = pk_idx;
while lo > 1 && pwr_norm(lo-1) >= 0.5
    lo = lo - 1;
end
hi = pk_idx;
while hi < N_dop && pwr_norm(hi+1) >= 0.5
    hi = hi + 1;
end

bw_3dB = abs(v_axis(hi) - v_axis(lo));
bw_3dB = max(bw_3dB, 0.5);      % floor 0.5 m/s
end


% ─────────────────────────────────────────────────────────────────────────
function [rd_out, diag] = tf_ridge_filter(frm, win_doppler, tfc, dv_att_pulse, sigma_pred_v, tgt_range_bin, B_flat_fix_ms)
%  慢时间STFT + Viterbi脊线追踪 + 目标保持型时频掩模 + WOLA-ISTFT
%  frm: [N_slow × N_fast] MOCO+平动补偿后的 慢时间×距离
%  tgt_range_bin: 目标距离门(仅用于导出目标距离区时频谱诊断, 不参与滤波)
%  返回 rd_out: 时频滤波后再做多普勒FFT 的 RD 图 [N_slow × N_fast]
[N_slow, N_fast] = size(frm);
L=tfc.L; H=tfc.H; Nfft=tfc.Nfft; pad=tfc.pad; w=tfc.w(:); v_tf=tfc.v_tf(:);
dv_bin = v_tf(2)-v_tf(1);

% ── 1) 慢时间 STFT（补零覆盖整 CPI，向量化跨距离门）──
xp = [zeros(pad,N_fast); frm; zeros(pad,N_fast)];
Np = size(xp,1);
starts = 1:H:(Np-L+1);
Nf = numel(starts);
X_tf = zeros(Nfft, Nf, N_fast);
Ptf  = zeros(Nfft, Nf);                       % 杂波距离区时频功率（脊线追踪用）
Ptf_tgt = zeros(Nfft, Nf);                    % 目标距离区时频功率（诊断/标注用）
if nargin>=6 && ~isempty(tgt_range_bin)
    tr = max(1,tgt_range_bin-1):min(N_fast,tgt_range_bin+1);
else
    tr = [];
end
for i = 1:Nf
    s = starts(i);
    seg = xp(s:s+L-1,:) .* w;                 % [L×N_fast]
    F   = fftshift(fft(seg, Nfft, 1), 1);     % [Nfft×N_fast]
    X_tf(:,i,:) = F;
    Ptf(:,i) = sum(abs(F(:,tfc.clut_bins)).^2, 2);
    if ~isempty(tr)
        Ptf_tgt(:,i) = sum(abs(F(:,tr)).^2, 2);
    end
end

% ── 2) 脊线预测：数据锚点 + IMU 蜿蜒先验 ──
Ptf_sum = sum(Ptf,2);
[~, k_anchor] = max(Ptf_sum);                 % 数据驱动锚定杂波中心 bin
% 把逐脉冲 IMU 径向速度映射到每个 STFT 帧，仅取去均值的"蜿蜒形状"
delta_v = zeros(Nf,1);
for n = 1:Nf
    c0 = starts(n) - pad + L/2;               % 帧中心在慢时间坐标
    a = min(max(round(c0-L/2),1),N_slow);
    b = min(max(round(c0+L/2),1),N_slow);
    if b<a; tmp=a; a=b; b=tmp; end
    delta_v(n) = mean(dv_att_pulse(a:b));
end
delta_v = delta_v - mean(delta_v);
k_pred  = round(k_anchor + delta_v/dv_bin);
k_pred  = min(max(k_pred,1),Nfft);

% ── 3) Viterbi 脊线追踪 ──
k_ridge = viterbi_ridge(Ptf, k_pred, tfc.dk_search, tfc.sigmaT, tfc.lamT);

% ── 4) 逐帧 3dB 谱宽（数据）与 IMU 预测谱宽融合 → 平顶/滚降半宽(bin)──
sigma_pred_bin = max(sigma_pred_v/dv_bin, 1);
B_flat_bin = zeros(Nf,1);
B_roll_bin = zeros(Nf,1);
for n = 1:Nf
    col = Ptf(:,n);
    pk = col(k_ridge(n))+eps;
    above = col >= 0.1*pk;
    bw = local_3db_width(above, k_ridge(n));
    sig_data_bin = max(bw/1.5, 2.0);
    sig_bin = sqrt(sig_data_bin^2 + sigma_pred_bin^2);
    B_flat_bin(n) = min(tfc.kappa*sig_bin, tfc.Bflat_cap/dv_bin);
    B_roll_bin(n) = 0.5*sig_bin;
end

% ── 5) 目标保持型时频掩模 G(k,n)：平顶 + 高斯滚降（H_min 兜底，不整行保护）──
%  注：目标保护改由主程序在输出 RD 上做 2D 局部恢复，避免整行抬增益造成零频亮线
    if nargin>=7 && ~isempty(B_flat_fix_ms) && B_flat_fix_ms>0
        B_flat_bin(:) = B_flat_fix_ms / dv_bin;   % 非自适应：固定平顶宽度
    end
G  = ones(Nfft, Nf);
kk = (1:Nfft).';
for n = 1:Nf
    d = abs(kk - k_ridge(n));
    g = ones(Nfft,1);
    flat = d <= B_flat_bin(n);
    g(flat) = tfc.Hmin;
    roll = ~flat;
    g(roll) = max(tfc.Hmin, 1 - tfc.alpha*exp(-((d(roll)-B_flat_bin(n)).^2)/(2*max(B_roll_bin(n),0.5)^2)));
    G(:,n) = g;
end

% ── 6a) 距离自适应 CNR 权重 w_r ──
%  思路：杂波能量集中在脊线附近的速度带内（in-band），带外基本是噪声。
%  用带内 excess power / 带外噪声估计 CNR(r)，转化为杂波存在权重 w(r)∈[0,1]。
%  w(r)≈1 时（强杂波）全幅度使用掩模；w(r)≈0 时（无杂波）掩模不生效, 避免对纯噪声陷波。
in_band_kv  = abs((1:Nfft).' - round(mean(k_ridge))) <= round(mean(B_flat_bin));  % [Nfft×1]
out_band_kv = ~in_band_kv;
n_in = sum(in_band_kv); n_out = max(sum(out_band_kv),1);
P_in_r  = zeros(1, N_fast);
P_out_r = zeros(1, N_fast);
for i = 1:Nf
    F_i = squeeze(X_tf(:,i,:));          % [Nfft×N_fast]
    P_in_r  = P_in_r  + sum(abs(F_i(in_band_kv,:)).^2, 1);
    P_out_r = P_out_r + sum(abs(F_i(out_band_kv,:)).^2, 1);
end
P_out_norm   = P_out_r * (n_in/n_out);  % 归一化到相同带宽
CNR_excess   = max(0, P_in_r - P_out_norm) ./ (P_out_norm + eps);
w_r          = CNR_excess ./ (1 + CNR_excess);   % [1×N_fast]∈[0,1]
w_r          = movmean(w_r, 7);                   % 距离平滑 7 点
w_r          = min(max(w_r, 0), 1);

% ── 6b) 掩模 + WOLA-ISTFT 重建（距离自适应掩模）──
yp = zeros(Np, N_fast); ws = zeros(Np,1);
for i = 1:Nf
    s = starts(i);
    F_i  = squeeze(X_tf(:,i,:));         % [Nfft×N_fast]
    % 距离自适应掩模: G_r(k,r) = 1 + w(r)*(G(k)-1)
    % 无杂波处(w≈0): G_r≈1, 不抑制; 强杂波处(w≈1): G_r=G, 全幅抑制
    G_r  = 1 + w_r .* (G(:,i) - 1);    % [Nfft×N_fast], broadcast
    col  = F_i .* G_r;
    seg  = ifft(ifftshift(col, 1), Nfft, 1);
    seg  = seg(1:L,:);
    yp(s:s+L-1,:) = yp(s:s+L-1,:) + seg .* w;
    ws(s:s+L-1)   = ws(s:s+L-1) + w.^2;
end
ws(ws<1e-12) = 1e-12;
frm_filt = yp(pad+1:pad+N_slow,:) ./ ws(pad+1:pad+N_slow);

% ── 7) 多普勒 FFT → RD ──
rd_out = fftshift(fft(frm_filt .* win_doppler, N_slow, 1), 1);

% ── 诊断（可视化用）──
diag.Ptf=Ptf; diag.Ptf_tgt=Ptf_tgt; diag.Ptf_masked=Ptf.*G.^2;
diag.k_pred=k_pred; diag.k_ridge=k_ridge; diag.G=G; diag.Nf=Nf;
diag.v_ridge = v_tf(k_ridge); diag.B_flat_v = B_flat_bin*dv_bin; diag.B_roll_v = B_roll_bin*dv_bin;
diag.w_r = w_r; diag.CNR_excess = CNR_excess;
% 【断裂点1 修正】暴露 TF 掩模实际使用的平顶半宽（转回 m/s）。
%   原代码主循环用 RD 域估计的 B_flat_final 做速度门控，而 TF 掩模内部独立重算了
%   一套更宽的 B_flat_bin，两者相差 2~3 倍且从不同步。这里把内部真正使用的宽度
%   导出，主循环改用它做速度门控，使"门控排除带"与"TF 掩模压制带"对齐。
diag.B_flat_used_ms = mean(B_flat_bin) * dv_bin;   % 平均平顶半宽 (m/s)
diag.B_flat_max_ms  = max(B_flat_bin)  * dv_bin;   % 最大平顶半宽 (m/s)，更保守
end

% ─────────────────────────────────────────────────────────────────────────
function ridge = viterbi_ridge(Ptf, k_pred, dk, sigmaT, lamT)
%VITERBI_RIDGE  IMU 先验约束 + 能量驱动的时频脊线动态规划追踪
[Kf, Nf] = size(Ptf);
kk = (1:Kf).';
sigp = max(dk/2, 1);
E = log(Ptf + 1e-12) - ((kk - k_pred(:).').^2)/(2*sigp^2);   % 数据项 + 高斯先验
win = abs(kk - k_pred(:).') <= dk;
NEG = -1e18;
E(~win) = NEG;
V = NEG*ones(Kf, Nf); bp = ones(Kf, Nf);
V(:,1) = E(:,1);
for n = 2:Nf
    Vprev = V(:,n-1);
    for k = 1:Kf
        if ~win(k,n); continue; end
        trans = -((k - kk).^2)/(2*sigmaT^2)*lamT;
        [val, kb] = max(Vprev + trans);
        V(k,n) = val + E(k,n);
        bp(k,n) = kb;
    end
end
ridge = zeros(Nf,1);
[~, ridge(Nf)] = max(V(:,Nf));
for n = Nf-1:-1:1
    ridge(n) = bp(ridge(n+1), n+1);
end
end

% ─────────────────────────────────────────────────────────────────────────
function bw = local_3db_width(above, kc)
%LOCAL_3DB_WIDTH  围绕脊线 kc 的连续 -3dB 带宽（bin）
N = numel(above);
if ~above(kc)
    bw = 1;
    return;
end
lo = kc;
while lo>1 && above(lo-1)
    lo=lo-1;
end
hi = kc;
while hi<N && above(hi+1)
    hi=hi+1;
end
bw = hi - lo + 1;
end

% ─────────────────────────────────────────────────────────────────────────
function H = adaptive_soft_notch(v_axis, mu_c, B_flat, B_roll, alpha, H_min, N_dop, N_rng)
%ADAPTIVE_SOFT_NOTCH  目标失真约束的自适应软陷波（平顶+高斯滚降）—— 沿用 v2
dv = abs(v_axis(:) - mu_c);
H_vec = ones(N_dop,1);
flat_mask = dv <= B_flat;  H_vec(flat_mask) = H_min;
roll_idx = find(~flat_mask);
if ~isempty(roll_idx) && B_roll>0
    dv_roll = dv(roll_idx) - B_flat;
    H_vec(roll_idx) = max(H_min, 1 - alpha*exp(-dv_roll.^2/(2*B_roll^2)));
end
H = repmat(H_vec, 1, N_rng);
end

% ─────────────────────────────────────────────────────────────────────────
function H = make_wide_notch(v_axis, center, half_bw, N_dop, N_rng)
%MAKE_WIDE_NOTCH  硬零陷滤波器
H_vec = ones(N_dop,1);
H_vec(abs(v_axis(:)-center) <= half_bw) = 0;
H = repmat(H_vec, 1, N_rng);
end

% ─────────────────────────────────────────────────────────────────────────
function R = build_R_zyx_local(phi, theta, psi)
%BUILD_R_ZYX_LOCAL  ZYX 欧拉角旋转矩阵
Rz = [cos(psi),-sin(psi),0; sin(psi),cos(psi),0; 0,0,1];
Ry = [cos(theta),0,sin(theta); 0,1,0; -sin(theta),0,cos(theta)];
Rx = [1,0,0; 0,cos(phi),-sin(phi); 0,sin(phi),cos(phi)];
R  = Rz*Ry*Rx;
end

% ─────────────────────────────────────────────────────────────────────────
function s_clut = sim_dynamic_clutter(Num_chirp, N_fast, ~, fc, PRI, K, Tc, ...
    t_fast, ant_pos_all, ant_vel_all, lambda, tgt_coor0, SCR_db)
%SIM_DYNAMIC_CLUTTER  物理散射单元模型的动态地面杂波（沿用 v2）
N_scatter = 80; c = 3e8; tgt_R0 = norm(tgt_coor0);
rng(42);
r_scatter  = linspace(50, 900, N_scatter);
az_scatter = deg2rad(linspace(-30, 30, N_scatter) + 5*randn(1, N_scatter));
sigma0     = randn(1, N_scatter) + 0.5j*randn(1, N_scatter);
s_clut = zeros(Num_chirp, N_fast);
for i = 1:Num_chirp
    t_curr  = (i-1)*PRI; ant_pos = ant_pos_all(i,:); ant_vel = ant_vel_all(i,:);
    for n = 1:N_scatter
        scat_pos = [r_scatter(n)*cos(az_scatter(n)), r_scatter(n)*sin(az_scatter(n)), 0];
        dist_sn  = norm(scat_pos - ant_pos);
        e_sn     = (scat_pos - ant_pos)/dist_sn;
        fd_sn    = -2/lambda*dot(ant_vel, e_sn);
        tau_sn   = 2*dist_sn/c;
        amp_sn   = (tgt_R0/dist_sn)^2 * sigma0(n);
        s_clut(i,:) = s_clut(i,:) + amp_sn*rectpuls(t_fast-Tc/2-tau_sn, Tc) .* ...
            exp(1j*pi*K*(t_fast-tau_sn).^2) .* exp(-1j*2*pi*fc*tau_sn) .* exp(1j*2*pi*fd_sn*t_curr);
    end
end
P_clut = mean(abs(s_clut).^2,'all'); scr_lin = 10^(SCR_db/10);
s_clut = s_clut * sqrt(1/(scr_lin*P_clut + eps));
end

% ─────────────────────────────────────────────────────────────────────────
function [mu_c, sigma_c] = estimate_clutter_stats_rd(rd, v_axis, mu_prior, search_hw, clut_bins)
%  从 RD 图估计残余杂波中心与谱宽（沿用 v2，断裂点2 修正）
%  【断裂点2 修正】原代码 pwr = mean(abs(rd).^2, 2) 对全部 N_fast 个距离门取均值。
%    真实杂波只分布在近距离的杂波距离区，大量远距离门只有噪声；均值被噪声稀释后
%    杂波谱"峰"被压低、-10 dB 带宽收缩，导致 σ_c 系统性偏小。
%    修正：仅对杂波主导的距离区 clut_bins 内的门做平均。
if nargin<3 || isempty(mu_prior);  mu_prior = 0;   end
if nargin<4 || isempty(search_hw); search_hw = 4.0; end
v_axis = v_axis(:);
if nargin>=5 && ~isempty(clut_bins)
    cb  = clut_bins(clut_bins>=1 & clut_bins<=size(rd,2));   % 防越界
    if isempty(cb); cb = 1:size(rd,2); end
    pwr = mean(abs(rd(:, cb)).^2, 2);   % 仅杂波距离区取均值
else
    pwr = mean(abs(rd).^2, 2);          % 兼容旧调用：全距离门
end
in_win = abs(v_axis - mu_prior) <= search_hw;
if nnz(in_win) < 5
    mu_c = mu_prior;
    sigma_c = 1.0;
    return;
end
v_w = v_axis(in_win);
p_w = pwr(in_win);
p_w = p_w - min(p_w);
pk = max(p_w)+eps;
band = p_w >= 0.1*pk;
v_m = v_w(band);
p_m = p_w(band);
p_sum = sum(p_m)+eps;
mu_c = sum(v_m.*p_m)/p_sum;
sigma_c = sqrt(sum((v_m-mu_c).^2.*p_m)/p_sum);
end

% ─────────────────────────────────────────────────────────────────────────
function model = compute_clutter_ridge_model(vel_all, ant_vel_all, ~, ...
    pos_all, ~, ~, ~, N_slow, step, num_frames, Num_chirp, radar_offset, ...
    beam_dir, theta_bw, IN_coor0, ~)
%COMPUTE_CLUTTER_RIDGE_MODEL  基于 IMU 运动状态的每 CPI 脊线参数（沿用 v2）
model.mu_c_plat   = zeros(num_frames,1);
model.mu_c_att    = zeros(num_frames,1);
model.sigma_c_est = zeros(num_frames,1);
model.clutter_pwr_ratio = zeros(num_frames,1);
for f = 1:num_frames
    s_idx = (f-1)*step + 1; e_idx = min(s_idx+N_slow-1, Num_chirp);
    frame_idx = s_idx:e_idx;
    v_mean = mean(vel_all(frame_idx,:),1); v_mag_plat = norm(v_mean);
    model.mu_c_plat(f) = v_mag_plat;
    mid_idx = min(max(round((s_idx+e_idx)/2),1), Num_chirp);
    apc_mid = pos_all(mid_idx,:) + radar_offset;
    los_vec = IN_coor0 - apc_mid; e_los = los_vec/norm(los_vec);
    v_att_all_frame = ant_vel_all(frame_idx,:) - vel_all(frame_idx,:);
    vr_att = v_att_all_frame * e_los.';
    model.mu_c_att(f) = mean(vr_att);
    sigma_c_cpi = std(vr_att);
    heading_f = atan2(v_mean(2), v_mean(1));
    squint_f  = beam_dir - heading_f;
    sigma_c_beam = v_mag_plat*abs(sin(squint_f))*theta_bw/(2*sqrt(3));
    % 注：sigma_icm 参数保留接口但不再加入，由主循环数据驱动估计后动态叠加
    % 【断裂点3 修正】车载前视雷达即使 squint≈0（波束与运动方向几乎平行），
    %   ICM（地面内禀运动：植被晃动、地表微动等）仍会造成 0.5~1.0 m/s 的展宽。
    %   原代码 squint≈0 时 sigma_c_beam 仅 ~0.1 m/s，导致先验严重低估、陷波过窄。
    %   这里加入保守下界 SIGMA_C_FLOOR，使先验不低于真实杂波的物理展宽量级。
    SIGMA_C_FLOOR = 0.8;   % 杂波先验展宽保守下界 (m/s)，覆盖 ICM 等几何先验未建模的展宽
    sigma_c = max(sqrt(sigma_c_beam^2 + sigma_c_cpi^2), SIGMA_C_FLOOR);
    model.sigma_c_est(f) = sigma_c;    % 平台几何分量 + ICM 保守下界
    model.clutter_pwr_ratio(f) = min(1.0, v_mag_plat/15);
end
end

% ─────────────────────────────────────────────────────────────────────────
function [d_found, r_found, d_refined, r_refined, found_ok] = find_target_peak_prior(rd, d_prior, r_prior, W_d, W_r, gamma_min_dB)
%FIND_TARGET_PEAK_PRIOR  先验引导的目标峰值搜索（沿用 v2）
[N_dop, N_rng] = size(rd);
d_lo=max(1,d_prior-W_d); d_hi=min(N_dop,d_prior+W_d);
r_lo=max(1,r_prior-W_r); r_hi=min(N_rng,r_prior+W_r);
win_pwr = abs(rd(d_lo:d_hi, r_lo:r_hi)).^2;
[~, pk] = max(win_pwr(:)); [dr, rr] = ind2sub(size(win_pwr), pk);
d_found = d_lo+dr-1; r_found = r_lo+rr-1;
P_peak = win_pwr(dr,rr); P_med = median(win_pwr(:))+eps;
if 10*log10(P_peak/P_med) >= gamma_min_dB
    found_ok = true;
else
    d_found=d_prior; r_found=r_prior; found_ok=false;
end
d_refined=double(d_found); r_refined=double(r_found);
end

% ─────────────────────────────────────────────────────────────────────────
function [scr_dB, P_tgt] = compute_scr(rd, tgt_dop, tgt_bin, N_dop, N_rng, clut_bins)
%COMPUTE_SCR  信杂比 (dB)
%  P_tgt ：以峰值定位位置为中心取 ±1 bin 窗口均值
%          （避免原 max(7×7) 可能捕获残余杂波尖峰的问题）
%  P_c   ：仅在杂波主导距离区（clut_bins）计算零多普勒 ±25 bin 平均功率
%          （避免远距离纯噪声门稀释杂波功率，使 SCR 虚高）
if nargin < 6 || isempty(clut_bins)
    clut_bins = 1:N_rng;   % 兼容旧调用：默认全距离门
end
rb = 1; db = 1;
r_idx = max(1,tgt_bin-rb):min(N_rng,tgt_bin+rb);
d_idx = max(1,tgt_dop-db):min(N_dop,tgt_dop+db);
P_tgt = mean(abs(rd(d_idx,r_idx)).^2, 'all');
dc = round(N_dop/2); cld = max(1,dc-25):min(N_dop,dc+25);
cb = clut_bins(clut_bins>=1 & clut_bins<=N_rng);
if isempty(cb); cb = 1:N_rng; end
P_c = mean(abs(rd(cld,cb)).^2, 'all');
scr_dB = 10*log10(P_tgt/(P_c+eps));
end

% ─────────────────────────────────────────────────────────────────────────
function bw = compute_bw3dB(rd, clut_bins, v_axis, N_dop)
%COMPUTE_BW3DB  杂波多普勒 -3dB 带宽（m/s）（沿用 v2）
pwr = zeros(N_dop,1);
for b = clut_bins
    if b <= size(rd,2); pwr = pwr + abs(rd(:,b)).^2; end
end
pwr_n = pwr/(max(pwr)+eps);
bw = sum(pwr_n>=0.1)*(v_axis(2)-v_axis(1));
end

% ─────────────────────────────────────────────────────────────────────────
function [det, Z_map] = my_go_cfar(pwr_map, Tr, Gr, Td, Gd, Pfa)
%MY_GO_CFAR  标准 GO-CFAR 检测（向量化，基于卷积）
%  修正：
%    - Z_d_avg = sum(2*Td 参考单元) / (2*Td)  多普勒维平均
%    - Z_r_avg = sum(2*Tr 参考单元) / (2*Tr)  距离维平均
%    - Z_map   = max(Z_d_avg, Z_r_avg)        GO = max of averages
%    - N_ref   = 2*(Td+Tr)                    总参考单元数
%    - alpha   = N_ref * (Pfa^(-1/N_ref)-1)
%    - det     = (pwr > alpha * Z_map)
[N_dop, N_rng] = size(pwr_map);
N1 = 2*Td; N2 = 2*Tr;            % 各方向参考单元数
N_ref = N1 + N2;                  % 总参考单元数 = 2*(Td+Tr)

% 多普勒维参考窗：上下各 Td 个，隔 Gd 保护
kernel_dop = [ones(Td,1); zeros(2*Gd+1,1); ones(Td,1)];
Z_d_avg = conv2(pwr_map, kernel_dop, 'same') / N1;

% 距离维参考窗：左右各 Tr 个，隔 Gr 保护
kernel_rng = [ones(1,Tr), zeros(1,2*Gr+1), ones(1,Tr)];
Z_r_avg = conv2(pwr_map, kernel_rng, 'same') / N2;

% GO = Greatest Of (max of averages)
Z_map = max(Z_d_avg, Z_r_avg);
Z_map(Z_map < eps) = eps;

% 有效区域
valid_mask = false(N_dop, N_rng);
dv_lo = Td+Gd+1; dv_hi = N_dop-Td-Gd;
rv_lo = Tr+Gr+1; rv_hi = N_rng-Tr-Gr;
valid_mask(dv_lo:dv_hi, rv_lo:rv_hi) = true;
Z_map(~valid_mask) = NaN;

% 标准 GO-CFAR 检测门限
alpha = N_ref * (Pfa^(-1/N_ref) - 1);
det = (pwr_map > alpha * Z_map) & valid_mask;
end

% ─────────────────────────────────────────────────────────────────────────
function [cand_mask, cand_list] = detect_candidate_targets(rd_power, thresh_db)
%DETECT_CANDIDATE_TARGETS  局部峰值检测 + 非极大值抑制
[N_dop, N_rng] = size(rd_power);
win = 10;
pad_power = padarray(rd_power, [win win], 'replicate', 'both');
local_mean = zeros(N_dop, N_rng);
for i = 1:N_dop
    for j = 1:N_rng
        window = pad_power(i:i+2*win, j:j+2*win);
        local_mean(i,j) = median(window(:));
    end
end
peak_ratio = rd_power ./ (local_mean + eps);
cand_mask = peak_ratio > 10^(thresh_db/10);
cand_mask = imregionalmax(rd_power) & cand_mask;
cand_list = find(cand_mask);
end

% ─────────────────────────────────────────────────────────────────────────
function is_isolated = isolated_peak_test(rd_power, r_idx, d_idx, guard_win, ref_win)
%ISOLATED_PEAK_TEST  孤立点目标判定（局部对比 + 距离延展性，剔除杂波脊）（沿用 v2）
[N_dop, N_rng] = size(rd_power);
r_lo_g=max(1,r_idx-guard_win); r_hi_g=min(N_rng,r_idx+guard_win);
d_lo_g=max(1,d_idx-guard_win); d_hi_g=min(N_dop,d_idx+guard_win);
guard_mask = false(N_dop,N_rng); guard_mask(d_lo_g:d_hi_g, r_lo_g:r_hi_g)=true;
guard_mask(d_idx,r_idx)=false; P_guard = mean(rd_power(guard_mask));
r_lo_r=max(1,r_idx-ref_win); r_hi_r=min(N_rng,r_idx+ref_win);
d_lo_r=max(1,d_idx-ref_win); d_hi_r=min(N_dop,d_idx+ref_win);
ref_mask = false(N_dop,N_rng); ref_mask(d_lo_r:d_hi_r, r_lo_r:r_hi_r)=true;
ref_mask(d_lo_g:d_hi_g, r_lo_g:r_hi_g)=false; P_ref = mean(rd_power(ref_mask));
contrast_ok = (P_guard/(P_ref+eps)) > 4.0;
far=20; pk = rd_power(d_idx,r_idx)+eps;
idx_lo = max(1,r_idx-far-3):max(1,r_idx-far);
idx_hi = min(N_rng,r_idx+far):min(N_rng,r_idx+far+3);
far_cells = [rd_power(d_idx,idx_lo), rd_power(d_idx,idx_hi)];
not_a_ridge = (median(far_cells)/pk) < 0.5;
is_isolated = contrast_ok && not_a_ridge;
end