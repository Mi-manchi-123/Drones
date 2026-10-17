%% =========================================================================
%  ob_moco_v4_ablation.m
%  单通道车载雷达：完整消融实验平台 — 7 方案对比 + ROC 评估
%
%  支持从 base workspace 导入参数：
%    TGT_VT_OPTIONS : 目标速度向量 [vx, vy, vz] (默认 [3.2, 2.4, 0])
%    SNR_db_input   : 信噪比标量 dB (默认 -15)
% =========================================================================
close all; clc;

% ── 实验配置（由 run_reviewer_experiments.m 注入） ──
if ~exist('EXPERIMENT_CFG', 'var') || ~isstruct(EXPERIMENT_CFG)
    EXPERIMENT_CFG = struct();
end
if ~isfield(EXPERIMENT_CFG, 'batch_mode')
    EXPERIMENT_CFG.batch_mode = false;
end
if ~isfield(EXPERIMENT_CFG, 'random_seed')
    EXPERIMENT_CFG.random_seed = 42;
end
if ~isfield(EXPERIMENT_CFG, 'save_figures')
    EXPERIMENT_CFG.save_figures = ~EXPERIMENT_CFG.batch_mode;
end
if EXPERIMENT_CFG.batch_mode
    set(0, 'DefaultFigureVisible', 'off');
else
    set(0, 'DefaultFigureVisible', 'on');
end
rng(EXPERIMENT_CFG.random_seed, 'twister');
batch_mode = EXPERIMENT_CFG.batch_mode;
if isfield(EXPERIMENT_CFG, 'clutter_range_bins') && ~isempty(EXPERIMENT_CFG.clutter_range_bins)
    CLUTTER_RANGE_BINS_CFG = EXPERIMENT_CFG.clutter_range_bins;
else
    CLUTTER_RANGE_BINS_CFG = [];
end

% ── 导入外部参数（来自 run_v4_experiments 的注入） ──
%  注意：不要使用 clear all，否则会清掉调用脚本的工作空间变量
if isfield(EXPERIMENT_CFG, 'target_vt') && ~isempty(EXPERIMENT_CFG.target_vt)
    TGT_VT_OPTIONS = EXPERIMENT_CFG.target_vt;
elseif evalin('base', 'exist(''TGT_VT_OPTIONS'', ''var'')')
    TGT_VT_OPTIONS = evalin('base', 'TGT_VT_OPTIONS');
else
    TGT_VT_OPTIONS = [3.2, 2.4, 0];
end
if isfield(EXPERIMENT_CFG, 'snr_db')
    SNR_db = EXPERIMENT_CFG.snr_db;
elseif evalin('base', 'exist(''SNR_db_input'', ''var'')')
    SNR_db = evalin('base', 'SNR_db_input');
else
    SNR_db = -15;
end

% 清除本脚本可能遗留的工作空间变量（仅清除非注入变量）
clear SCR_data ROC_data peak_ok mu_data sigma_data;  

fprintf('============================================================\n');
fprintf('  完整消融实验平台 — 7 方案对比 + ROC 评估 (v4)\n');
fprintf('  SNR = %.0f dB | 目标速度 = (%.1f, %.1f, %.1f) m/s\n', ...
    SNR_db, TGT_VT_OPTIONS(1), TGT_VT_OPTIONS(2), TGT_VT_OPTIONS(3));
fprintf('============================================================\n\n');

%% =========================================================================
% PART 0：公共雷达参数（与 v2/v3 保持一致）
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
overlap  = 100;             % 帧间重叠
Num_chirp = floor(Ob_time / PRI);   % 总脉冲数 = 10000

if isempty(CLUTTER_RANGE_BINS_CFG)
    CLUTTER_RANGE_BINS = [50, min(N_fast, round(0.6*N_fast))];
else
    CLUTTER_RANGE_BINS = CLUTTER_RANGE_BINS_CFG;
end

% 距离轴
t_fast      = (0:N_fast-1) / fs;
range_axis  = t_fast * c / 2;
range_res   = c / (fs * 2);

% 速度轴（CPI 级 / FFT 后 fftshift 约定）—— 七方案统一使用 v_tf
step       = N_slow - overlap;
num_frames = floor((Num_chirp - N_slow) / step) + 1;
if isfield(EXPERIMENT_CFG, 'max_frames') && ~isempty(EXPERIMENT_CFG.max_frames)
    use_frames = max(1, min(num_frames, round(EXPERIMENT_CFG.max_frames)));
    num_frames = use_frames;
    Num_chirp = (num_frames - 1) * step + N_slow;
end
v_tf       = ((0:N_slow-1) - floor(N_slow/2)) * (1/PRI) / N_slow * lambda/2;   % m/s
v_tf       = v_tf(:).';

% 窗函数
win_range   = hanning(N_fast).';
win_doppler = hanning(N_slow);

fprintf('[参数] λ=%.4f m | 距离分辨率=%.2f m | 多普勒速度分辨率=%.3f m/s\n', ...
    lambda, range_res, v_tf(2)-v_tf(1));
fprintf('[参数] 总脉冲数=%d | CPI帧数=%d | 每帧 N_slow=%d\n\n', ...
    Num_chirp, num_frames, N_slow);

%% =========================================================================
% PART 1：平台运动状态生成（ESKF数据 或 仿真替代）
% =========================================================================
USE_SIM_IMU = true;

if isfield(EXPERIMENT_CFG, 'imu_file') && ~isempty(EXPERIMENT_CFG.imu_file) && ...
        exist(EXPERIMENT_CFG.imu_file, 'file')
    fprintf('[数据] 加载实验指定 IMU 文件: %s\n', EXPERIMENT_CFG.imu_file);
    eskf_data = readtable(EXPERIMENT_CFG.imu_file);
    USE_SIM_IMU = false;
elseif exist('Data_deduced_in_reverse.csv', 'file')
    fprintf('[数据] 加载 ESKF 反推文件...\n');
    eskf_data = readtable('Data_deduced_in_reverse.csv');
    USE_SIM_IMU = false;
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
% PART 2：目标参数
% =========================================================================
fprintf('[目标] 配置目标参数...\n');
if isfield(EXPERIMENT_CFG, 'multi_targets') && ~isempty(EXPERIMENT_CFG.multi_targets)
    tgt = EXPERIMENT_CFG.multi_targets;
else
    tgt(1).coor0 = [300, 200, 100];
    tgt(1).vt    = TGT_VT_OPTIONS;
    tgt(1).A     = [0, 0, 0];
    tgt(1).type  = 'const';
end
N_tgt = numel(tgt);
if isfield(EXPERIMENT_CFG, 'target_index') && ~isempty(EXPERIMENT_CFG.target_index)
    target_index = min(N_tgt, max(1, round(EXPERIMENT_CFG.target_index)));
else
    target_index = 1;
end
tgt_sel = tgt(target_index);
if isfield(EXPERIMENT_CFG, 'radar_offset') && ~isempty(EXPERIMENT_CFG.radar_offset)
    radar_offset = EXPERIMENT_CFG.radar_offset;
else
    radar_offset = [0, 0, 5];   % 雷达阵面相对车体质心的安装偏移 (m)
end
fprintf('       目标初始位置=(%.0f,%.0f,%.0f)m  速度=(%.1f,%.1f,%.1f)m/s\n', ...
    tgt_sel.coor0(1), tgt_sel.coor0(2), tgt_sel.coor0(3), ...
    tgt_sel.vt(1), tgt_sel.vt(2), tgt_sel.vt(3));

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
% PART 4：生成目标回波
% =========================================================================
IN_coor0 = mean(reshape([tgt.coor0], 3, []), 2)';
tgt_R0   = norm(IN_coor0 - ant_pos_all(1,:));

s_rx_raw    = zeros(Num_chirp, N_fast);
target_pos  = zeros(Num_chirp, 3, N_tgt);
dist_target_arr = zeros(Num_chirp, N_tgt);

for i = 1:Num_chirp
    t_curr       = (i-1) * PRI;
    ant_pos      = ant_pos_all(i,:);
    for k = 1:N_tgt
        tgt_pos_k = tgt(k).coor0 + t_curr * tgt(k).vt + 0.5 * tgt(k).A * t_curr^2;
        target_pos(i,:,k) = tgt_pos_k;
        dist_t = norm(tgt_pos_k - ant_pos);
        dist_target_arr(i,k) = dist_t;
        range_amp = (tgt_R0 / dist_t)^2;
        tau = 2 * dist_t / c;
        s_rx_raw(i,:) = s_rx_raw(i,:) + ...
            range_amp * rectpuls(t_fast - Tc/2 - tau, Tc) .* ...
            exp(1j*pi*K*(t_fast - tau).^2) .* exp(-1j*2*pi*fc*tau);
    end
end
target_pos_sel = target_pos(:,:,target_index);

%% =========================================================================
% PART 5：生成或加载杂波信号
% =========================================================================
USE_SIM_CLUTTER = true;

if isfield(EXPERIMENT_CFG, 'clutter_file') && ~isempty(EXPERIMENT_CFG.clutter_file) && ...
        exist(EXPERIMENT_CFG.clutter_file, 'file')
    fprintf('[数据] 加载实验指定杂波数据: %s\n', EXPERIMENT_CFG.clutter_file);
    data = load(EXPERIMENT_CFG.clutter_file, 's_clutter_full');
    s_clutter_full = data.s_clutter_full;
    USE_SIM_CLUTTER = false;
elseif exist('s_clutter_full.mat', 'file')
    fprintf('[数据] 加载预生成杂波数据...\n');
    data = load('s_clutter_full.mat', 's_clutter_full');
    s_clutter_full = data.s_clutter_full;
    USE_SIM_CLUTTER = false;
else
    fprintf('[仿真] 生成动态地面杂波（散射单元模型）...\n');
    s_clutter_full = sim_dynamic_clutter(Num_chirp, N_fast, fs, fc, PRI, K, Tc, ...
        t_fast, ant_pos_all, ant_vel_all, lambda, ...
        IN_coor0, -10);
    % 缓存杂波到文件，后续 SNR/速度组合直接加载（杂波不依赖 SNR 和目标速度）
    save('s_clutter_full.mat', 's_clutter_full', '-v7.3');
    fprintf('[仿真] 杂波生成完毕，已缓存至 s_clutter_full.mat。\n\n');
end
if size(s_clutter_full, 1) > Num_chirp
    s_clutter_full = s_clutter_full(1:Num_chirp, :);
end

%% =========================================================================
% PART 6：运动补偿（MOCO）+ 脉冲压缩
% =========================================================================
dist1 = vecnorm(ant_pos_all  - target_pos_sel, 2, 2);    % 含姿态扰动瞬时斜距
dist2 = vecnorm(ant_ideal_all - target_pos_sel, 2, 2);   % 理想瞬时斜距
delta_swing_all     = dist1 - dist2;
swing_phase_2d_corr = exp(1j * 4*pi * delta_swing_all / lambda);

t_ref = (0:N_fast-1) / fs;
s_t   = rectpuls(t_ref - Tc/2, Tc) .* exp(1j*pi*K*t_ref.^2);
X_tx  = fft(s_t, N_fast);

% ── 关键修正：SNR 作为外部参数，主循环外完成信号生成 ──
%   SNR_db 将由外层实验脚本循环传入；此处使用默认值 -15 dB
%   （实际运行时可在脚本中修改）
% SNR_db 已在文件顶部导入 (默认 -15 dB)
% SNR_db       = -15;
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
% PART 7：创新点1 — 运动状态驱动的杂波脊线参数（每 CPI 标量预测）
% =========================================================================
fprintf('[创新1] 运动状态驱动的杂波脊线参数建模...\n');
ant0_approx = pos_all(1,:) + radar_offset;
beam_dir    = atan2(IN_coor0(2)-ant0_approx(2), IN_coor0(1)-ant0_approx(1));
theta_bw    = deg2rad(30);
clutter_model = compute_clutter_ridge_model(vel_all, ant_vel_all, ang_rate_all, ...
    pos_all, ang_all, lambda, PRI, N_slow, step, num_frames, Num_chirp, radar_offset, ...
    beam_dir, theta_bw, IN_coor0, 0);
fprintf('[创新1] 已计算 %d 帧平台几何杂波分量\n\n', num_frames);

%% =========================================================================
% PART 8：处理配置
% =========================================================================

% ── 通用融合参数（所有自适应方案共用）──
BETA_FUSION = 0.6;     % IMU 先验权重
KAPPA_BW    = 2.5;     % 陷波/掩模带宽扩展因子
ALPHA_EMA   = 0.7;     % EMA 帧间平滑
sigma_c_ema = 1.5;     % 运行时 EMA 初始状态
sigma_c_total_arr = zeros(num_frames,1);

% ── 方案A' 参数 ──
FIX_NOTCH_BW = 0.5;    % 固定零陷半宽 (m/s)

% ── 方案E (PCA) 参数 ──
PCA_N_COMPONENTS = 3;  % 杂波子空间秩（取前几个特征值）

% ── 时频处理配置（方案C/F/D/D-noIMU 共用）──
tfc.L         = 64;            % STFT 窗长
tfc.H         = 16;            % 帧移
tfc.Nfft      = N_slow;        % FFT 点数
tfc.pad       = tfc.L/2;       % 两端补零
tfc.w         = hann(tfc.L);   % Hann 窗
tfc.v_tf      = v_tf;
tfc.dk_search = 30;            % Viterbi 搜索半窗（bin）
tfc.sigmaT    = 3;             % 帧间转移平滑标准差（bin）
tfc.lamT      = 0.4;           % 转移项权重
tfc.kappa     = KAPPA_BW;
tfc.Hmin      = 0.01;          % 目标最低通过增益
tfc.alpha     = 1.0;           % 陷波中心最大深度
tfc.Bflat_cap = 4;             % 平顶半宽上限 (m/s)
BETA_CONSERVATIVE = 0.8;
DELTA_MARGIN = 0.1;

% ── 方案F (Hough) 参数 ──
houghOpts.n_theta   = 180;    % Hough 变换角度离散数
houghOpts.n_rho     = 200;    % Hough 变换距离离散数
houghOpts.npeaks    = 1;      % 检测脊线数
houghOpts.energy_th = 0.3;    % 二值化能量门限（相对峰值）

% ── CFAR 参数（ROC 评估用）──
Tr_go=8; Gr_go=2; Td_go=4; Gd_go=2;

fprintf('[配置] 融合 β=%.1f, κ=%.1f | STFT: L=%d, H=%d\n', ...
    BETA_FUSION, KAPPA_BW, tfc.L, tfc.H);
fprintf('       PCA 秩=%d | Hough θ=%d ρ=%d | 固定陷波 BW=%.1f m/s\n\n', ...
    PCA_N_COMPONENTS, houghOpts.n_theta, houghOpts.n_rho, FIX_NOTCH_BW);

%% =========================================================================
% PART 9：主处理循环 — 7 方案同步处理
% =========================================================================
fprintf('[主循环] 开始分帧处理（共 %d 帧，7 方案）...\n', num_frames);

% ── 方案记录结构体（7 方案 × 各指标）──
schemes = struct();
scheme_names = {'Ap', 'Bp', 'C', 'E', 'F', 'D', 'DnoIMU'};
n_schemes = 7;
for si = 1:n_schemes
    sn = scheme_names{si};
    schemes.(sn).SCR   = zeros(num_frames,1);
    schemes.(sn).TERI  = zeros(num_frames,1);
    schemes.(sn).FAR   = zeros(num_frames,1);
    schemes.(sn).Pd    = zeros(num_frames,1);
    schemes.(sn).BW3   = zeros(num_frames,1);
    schemes.(sn).pts   = zeros(num_frames,1);
end

% ── 公共记录 ──
peak_found  = false(num_frames,1);
mu_c_all    = zeros(num_frames,1);
sigma_c_all = zeros(num_frames,1);

% ── ROC 专用：预设 Pfa 列表 → 逐方案计算 ROC 曲线 ──
if isfield(EXPERIMENT_CFG, 'pfa_list') && ~isempty(EXPERIMENT_CFG.pfa_list)
    PFA_LIST = EXPERIMENT_CFG.pfa_list;
else
    PFA_LIST = logspace(-5, -1, 20);
end
roc_results = struct();
for si = 1:n_schemes
    sn = scheme_names{si};
    roc_results.(sn).Pd_vs_Pfa = zeros(num_frames, length(PFA_LIST));
    roc_results.(sn).FAR_vs_Pfa = zeros(num_frames, length(PFA_LIST));
end

% 中间帧保存（可视化用）
if isfield(EXPERIMENT_CFG, 'mid_frame') && ~isempty(EXPERIMENT_CFG.mid_frame)
    f_mid = min(num_frames, max(1, round(EXPERIMENT_CFG.mid_frame)));
else
    f_mid = min(30, num_frames);   % 帧数不足时使用可用的中间帧
end
MID = struct();

% ── 先验引导目标定位参数 ──
W_d_search = 4; W_r_search = 5; gamma_min_dB = 4.0;

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

    % ── 平台平动多普勒补偿（IMU）──
    e_los_bs   = R_fc;
    v_rad_plat = mean(vel_all(frame_idx,:) * e_los_bs, 1);
    fd_shift   = 2 * norm(v_rad_plat) / lambda;
    shift_plat = exp(-1j * 2*pi * fd_shift * slow_idx * PRI);

    frm_raw  = s_range_raw(frame_idx,:);
    frm_moco = s_range_moco(frame_idx,:);
    frm_moco_comp = frm_moco .* shift_plat;          % 平动+姿态MOCO（A'/B'/C/E/F/D/DnoIMU 共用）

    % ── 共用：MOCO 后未陷波 RD（残余杂波集中在 ~0）──
    rd_pre = fftshift(fft(frm_moco_comp .* win_doppler, N_slow, 1), 1);

    % ── 杂波中心/谱宽融合估计（创新点1+2）──
    mu_c_att_prior = clutter_model.mu_c_att(f_idx);
    clut_lo_est = max(1, min(N_fast, CLUTTER_RANGE_BINS(1)));
    clut_hi_est = max(clut_lo_est + 1, min(N_fast, CLUTTER_RANGE_BINS(2)));
    clut_bins_est = clut_lo_est:clut_hi_est;
    [mu_c_data, sigma_c_data] = estimate_clutter_stats_rd(rd_pre, v_tf, mu_c_att_prior, 4.0, clut_bins_est);
    mu_c_imu    = mu_c_att_prior;
    sigma_c_imu = clutter_model.sigma_c_est(f_idx);
    mu_c_final    = BETA_FUSION*mu_c_imu + (1-BETA_FUSION)*mu_c_data;
    mu_c_final    = max(min(mu_c_final, 1.5), -1.5);
    sigma_c_final = max(sigma_c_imu + sigma_c_data, 0.3);
    B_flat_final  = min(KAPPA_BW*sigma_c_final, tfc.Bflat_cap);
    B_roll_final  = 0.5*sigma_c_final;
    mu_c_all(f_idx) = mu_c_final;
    sigma_c_all(f_idx) = sigma_c_final;

    % ── 数据驱动总展宽估计 ──
    sigma_c_platform = clutter_model.sigma_c_est(f_idx);
    sigma_c_meas     = max(sigma_c_platform, sigma_c_data);
    sigma_c_ema      = ALPHA_EMA*sigma_c_ema + (1-ALPHA_EMA)*sigma_c_meas;
    sigma_c_ema      = max(sigma_c_ema, 0.8 * sigma_c_platform);
    sigma_c_total_arr(f_idx) = sigma_c_ema;
    sigma_c_final  = max(sigma_c_ema, 0.3);
    sigma_c_all(f_idx) = sigma_c_final;

    % ── 目标真值位置（仅评估用，不参与信号处理）──
    tp_n    = tgt_sel.coor0 + t_fc*tgt_sel.vt;
    los_vec = tp_n - ant_pos_fc;
    e_los_n = (los_vec/norm(los_vec)).';
    tgt_R_n = norm(los_vec);
    tgt_bin = min(max(round(tgt_R_n/range_res),1), N_fast);
    tgt_vr  = dot(tgt_sel.vt, e_los_n.');
    [~, tgt_dop_guess] = min(abs(v_tf - (-tgt_vr)));

    % ── 先验引导的目标峰值定位（所有方案共用同一参考位置）──
    [d_ref, r_ref, ~, ~, ok] = find_target_peak_prior(rd_pre, tgt_dop_guess, tgt_bin, ...
                                                       W_d_search, W_r_search, gamma_min_dB);
    peak_found(f_idx) = ok;
    if ~ok
        d_ref = tgt_dop_guess;
        r_ref = tgt_bin;
    end
    r_found = r_ref; d_found = d_ref;  % 统一别名

    % ── 杂波距离区（评估用）──
    clut_lo = clut_lo_est;
    clut_hi = max(clut_lo+1, tgt_bin-30);
    tfc.clut_bins = clut_lo:min(clut_hi, N_fast);
    clut_bins = tfc.clut_bins;

    % ════════════════════════════════════════════════════════════════
    %  方案 A' ：MOCO + 固定零陷 MTI（修正基线，隔离 MOCO 贡献）
    % ════════════════════════════════════════════════════════════════
    frm_Ap = frm_moco_comp .* win_doppler;   % 使用 MOCO 数据（与原 A 的关键区别）
    rd_Ap_pre = fftshift(fft(frm_Ap, N_slow, 1), 1);
    H_Ap  = make_wide_notch(v_tf, 0, FIX_NOTCH_BW, N_slow, N_fast);
    rd_Ap = rd_Ap_pre .* H_Ap;

    % ════════════════════════════════════════════════════════════════
    %  方案 B' ：MOCO + 自适应宽硬陷波
    %   σ_c 驱动带宽，而非固定 3.0 m/s
    % ════════════════════════════════════════════════════════════════
    B_hard_adaptive = KAPPA_BW * sigma_c_final;
    B_hard_adaptive = min(max(B_hard_adaptive, 0.5), 4.0);  % 限制范围
    H_Bp = make_wide_notch(v_tf, mu_c_final, B_hard_adaptive, N_slow, N_fast);
    rd_Bp = rd_pre .* H_Bp;

    % ════════════════════════════════════════════════════════════════
    %  方案 C ：MOCO + RD 域自适应软陷波（原 v3 方案 C）
    % ════════════════════════════════════════════════════════════════
    H_C  = adaptive_soft_notch(v_tf, mu_c_final, B_flat_final, B_roll_final, ...
                               tfc.alpha, tfc.Hmin, N_slow, N_fast);
    rd_C0 = rd_pre .* H_C;
    % 候选目标检测 → 2D 局部保护（与 v3 一致）
    cand_mask = detect_candidate_targets(abs(rd_C0).^2, 12);
    protect_cells = zeros(0,2);
    for k = find(cand_mask).'
        [d_idx, r_idx] = ind2sub([N_slow, N_fast], k);
        in_notch = abs(v_tf(d_idx) - mu_c_final) <= B_flat_final;
        if ~in_notch && isolated_peak_test(abs(rd_C0).^2, r_idx, d_idx, 5, 10)
            protect_cells(end+1,:) = [d_idx, r_idx];
        end
    end
    H_C2 = H_C;
    for ci = 1:size(protect_cells,1)
        dd = protect_cells(ci,1); rr = protect_cells(ci,2);
        H_C2(max(1,dd-2):min(N_slow,dd+2), max(1,rr-3):min(N_fast,rr+3)) = 1.0;
    end
    rd_C = rd_pre .* H_C2;

    % ════════════════════════════════════════════════════════════════
    %  方案 E ：MOCO + PCA 子空间投影杂波抑制
    % ════════════════════════════════════════════════════════════════
    rd_E = pca_clutter_suppression(frm_moco_comp, win_doppler, PCA_N_COMPONENTS, N_slow,N_fast);

    % ════════════════════════════════════════════════════════════════
    %  方案 F ：MOCO + STFT + Hough 变换脊线追踪 + 时频掩模
    % ════════════════════════════════════════════════════════════════
    % 使用与 D 方案相同的保守预压制 + Hough 脊线
    B_core_f = max(BETA_CONSERVATIVE * B_flat_final, 0.5);
    rd_pass1 = apply_conservative_notch(rd_pre, v_tf, mu_c_final, B_core_f);
    [rd_F, diagF] = hough_ridge_filter(frm_moco_comp, win_doppler, tfc, tgt_bin, houghOpts);

    % ════════════════════════════════════════════════════════════════
    %  方案 D ：MOCO + STFT + Viterbi 脊线追踪 + 时频掩模（本文）
    % ════════════════════════════════════════════════════════════════
    v_mean_f     = mean(vel_all(frame_idx,:),1);
    e_bore       = v_mean_f/(norm(v_mean_f)).';
    v_arm_f      = ant_vel_all(frame_idx,:) - vel_all(frame_idx,:);
    dv_att_pulse = -(v_arm_f .* e_bore);             % [N_slow×1]

    [rd_D_pure, diagD] = tf_ridge_filter(frm_moco_comp, win_doppler, tfc, ...
                                          dv_att_pulse, sigma_c_imu, tgt_bin);
    % D 方案候选目标检测 + 保护（同 v3）
    cand_mask_d = detect_candidate_targets(abs(rd_pass1).^2,10);
    protect_cells_D = zeros(0,2);
    n_cand_raw_D  = 0; n_cand_kept_D = 0;
    B_flat_for_gate = max(B_flat_final, diagD.B_flat_used_ms);
    for k = find(cand_mask_d).'
        [d_idx, r_idx] = ind2sub([N_slow, N_fast], k);
        n_cand_raw_D = n_cand_raw_D + 1;
        v_cand = v_tf(d_idx);
        in_clutter_zone = abs(v_cand - mu_c_final) <= (B_flat_for_gate + DELTA_MARGIN);
        if ~in_clutter_zone && isolated_peak_test(abs(rd_pass1).^2, r_idx, d_idx, 4, 8)
            protect_cells_D(end+1,:) = [d_idx, r_idx];
            n_cand_kept_D = n_cand_kept_D + 1;
        end
    end
    rd_D = rd_D_pure;
    for ci = 1:size(protect_cells_D,1)
        dd = protect_cells_D(ci,1); rr = protect_cells_D(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2);
        rrg = max(1,rr-3):min(N_fast,rr+3);
        rd_D(dr, rrg) = rd_pre(dr, rrg);
    end

    % ════════════════════════════════════════════════════════════════
    %  方案 D-noIMU ：Viterbi 中去除 IMU 先验
    % ════════════════════════════════════════════════════════════════
    [rd_Dno_pure, diagDno] = tf_ridge_filter_noimu(frm_moco_comp, win_doppler, tfc, ...
                                                    sigma_c_imu, tgt_bin);
    % D-noIMU 候选检测 + 保护（与 D 方案相同逻辑）
    cand_mask_dno = detect_candidate_targets(abs(rd_pass1).^2, 10);
    protect_cells_Dno = zeros(0,2);
    n_cand_raw_Dno  = 0; n_cand_kept_Dno = 0;
    B_flat_for_gate_dno = max(B_flat_final, diagDno.B_flat_used_ms);
    for k = find(cand_mask_dno).'
        [d_idx, r_idx] = ind2sub([N_slow, N_fast], k);
        n_cand_raw_Dno = n_cand_raw_Dno + 1;
        v_cand = v_tf(d_idx);
        in_clutter_zone = abs(v_cand - mu_c_final) <= (B_flat_for_gate_dno + DELTA_MARGIN);
        if ~in_clutter_zone && isolated_peak_test(abs(rd_pass1).^2, r_idx, d_idx, 4, 8)
            protect_cells_Dno(end+1,:) = [d_idx, r_idx];
            n_cand_kept_Dno = n_cand_kept_Dno + 1;
        end
    end
    rd_Dno = rd_Dno_pure;
    for ci = 1:size(protect_cells_Dno,1)
        dd = protect_cells_Dno(ci,1); rr = protect_cells_Dno(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2);
        rrg = max(1,rr-3):min(N_fast,rr+3);
        rd_Dno(dr, rrg) = rd_pre(dr, rrg);
    end

    % ════════════════════════════════════════════════════════════════
    %  性能指标统一计算（7 方案）
    % ════════════════════════════════════════════════════════════════
    rd_list = {rd_Ap, rd_Bp, rd_C, rd_E, rd_F, rd_D, rd_Dno};
    P_tgt_pre_ref = abs(rd_pre(d_found, r_found))^2 + eps;

    for si = 1:n_schemes
        sn = scheme_names{si};
        rd_s = rd_list{si};

        % SCR
        [schemes.(sn).SCR(f_idx), P_tgt_s] = compute_scr(rd_s, d_found, r_found, N_slow, N_fast);

        % TERI（统一分母：滤波前 rd_pre 的同一位置）
        schemes.(sn).TERI(f_idx) = P_tgt_s / P_tgt_pre_ref;

        % 杂波谱宽
        schemes.(sn).BW3(f_idx) = compute_bw3dB(rd_s, clut_bins, v_tf, N_slow);

        % ── CFAR 检测 + ROC 多 Pfa 统计（批处理：单次遍历算全部 21 个 Pfa）──
        pwr_s = abs(rd_s).^2;
        [det_map, totals_pfa, pd_pfa] = batch_go_cfar_roc(pwr_s, ...
            Tr_go, Gr_go, Td_go, Gd_go, 1e-4, PFA_LIST, d_found, r_found);
        schemes.(sn).pts(f_idx) = sum(det_map, 'all');
        schemes.(sn).Pd(f_idx)  = det_map(d_found, r_found);
        for pi = 1:length(PFA_LIST)
            roc_results.(sn).FAR_vs_Pfa(f_idx, pi) = totals_pfa(pi) / (N_slow * N_fast);
            roc_results.(sn).Pd_vs_Pfa(f_idx, pi)  = pd_pfa(pi);
        end
    end

    % ── 保存中间帧（可视化用）──
    if f_idx == f_mid
        MID.rd_Ap=rd_Ap; MID.rd_Bp=rd_Bp; MID.rd_C=rd_C;
        MID.rd_E=rd_E; MID.rd_F=rd_F;
        MID.rd_D=rd_D; MID.rd_Dno=rd_Dno; MID.rd_pre=rd_pre;
        MID.rd_pass1=rd_pass1; MID.mu_c=mu_c_final;
        MID.B_flat=B_flat_final; MID.B_roll=B_roll_final;
        MID.B_core=B_core_f; MID.tgt_bin=tgt_bin;
        MID.t_fc=t_fc; MID.v_tgt=-tgt_vr; MID.tgt_R=tgt_R_n;
        MID.d_ref=d_ref; MID.r_ref=r_ref; MID.ok=ok;
        MID.sigma_c_ema=sigma_c_ema;
        MID.clut_bins=clut_bins;
        MID.protect_cells_D=protect_cells_D;
        MID.n_cand_raw_D=n_cand_raw_D; MID.n_cand_kept_D=n_cand_kept_D;
        MID.B_flat_gate=B_flat_for_gate;
        MID.B_hard_adaptive=B_hard_adaptive;
        MID.diagD=diagD; MID.diagF=diagF; MID.diagDno=diagDno;
    end

    if mod(f_idx,20)==0
        fprintf('  帧 %3d/%d | μ_c=%.2f | SCR[C/D/Dno]=%.1f/%.1f/%.1f dB | 定位=%d\n', ...
            f_idx, num_frames, mu_c_final, ...
            schemes.C.SCR(f_idx), schemes.D.SCR(f_idx), schemes.DnoIMU.SCR(f_idx), ok);
    end
end
fprintf('\n[主循环] 处理完成。\n\n');

%% =========================================================================
% PART 10：性能指标汇总（7 方案消融对比）
% =========================================================================
N_cell = N_slow * N_fast;
vmask = peak_found; nv = sum(vmask);

fprintf('=============================================================\n');
fprintf('         性能指标汇总（%d 帧，有效定位帧 %d）\n', num_frames, nv);
fprintf('=============================================================\n');
fprintf('%-16s %9s %9s %9s %9s %9s %9s %9s\n', '指标', ...
    'A-prime', 'B-prime', 'C-RD', 'E-PCA', 'F-Hough', 'D-Viterbi', 'D-noIMU');
fprintf('%s\n', repmat('-',1,90));

for si = 1:n_schemes
    sn = scheme_names{si};
    scr_avg(si) = mean(schemes.(sn).SCR, 'omitnan');
    if nv>0
        teri_avg(si) = mean(schemes.(sn).TERI(vmask), 'omitnan');
    else
        teri_avg(si) = mean(schemes.(sn).TERI, 'omitnan');
    end
    bw3_avg(si) = mean(schemes.(sn).BW3(schemes.(sn).BW3>0), 'omitnan');
    far_avg(si) = mean(schemes.(sn).pts, 'omitnan')/N_cell;
end

fprintf('%-16s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n', '平均SCR(dB)', scr_avg);
if nv>0
    fprintf('%-16s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', '有效帧TERI', teri_avg);
end
fprintf('%-16s %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n', 'FAR(Pfa=1e-4)', far_avg);
fprintf('%-16s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n', 'SCR增益vsA-prime(dB)', ...
    scr_avg-scr_avg(1));
fprintf('%s\n', repmat('=',1,90));

% ── 消融链关键结论 ──
fprintf('\n[消融链分析]\n');
fprintf('  RD域链: A-prime -> B-prime -> C -> E\n');
fprintf('    SCR增益链: %.2f -> %.2f -> %.2f -> %.2f dB\n', scr_avg(1:4));
fprintf('  时频域链: F(Hough) -> D-noIMU -> D(Viterbi+IMU)\n');
fprintf('    SCR增益链: %.2f → %.2f → %.2f dB\n', scr_avg(5:7));
fprintf('\n  关键交叉对比:\n');
fprintf('    (1) 时频面 vs RD面: C(%.2f dB) vs D-noIMU(%.2f dB)  Δ=%.2f dB\n', ...
    scr_avg(3), scr_avg(7), scr_avg(7)-scr_avg(3));
fprintf('    (2) Viterbi vs Hough: F(%.2f dB) vs D-noIMU(%.2f dB)  Δ=%.2f dB\n', ...
    scr_avg(5), scr_avg(7), scr_avg(7)-scr_avg(5));
fprintf('    (3) IMU先验有用吗: D-noIMU(%.2f dB) vs D(%.2f dB)  Δ=%.2f dB\n', ...
    scr_avg(7), scr_avg(6), scr_avg(6)-scr_avg(7));
fprintf('\n');

%% =========================================================================
% PART 11：可视化（消融实验核心图）
% =========================================================================
if ~batch_mode
t_frames = (1:num_frames)*step*PRI;

% ── 图1：七方案中间帧 RD 对比 ──────────────────────────────────────────────
clim_rd = [-10, 100];
figure('Name','七方案RD对比','Position',[30,20,1580,760]);
rd_list_viz = {MID.rd_Ap, MID.rd_Bp, MID.rd_C, MID.rd_E, MID.rd_F, MID.rd_D, MID.rd_Dno};
ttl_viz = {'A-prime: MOCO+固定零陷(0.5m/s)', ...
           'B-prime: MOCO+自适应硬陷波', ...
           'C: MOCO+RD自适应软陷波', ...
           'E: MOCO+PCA子空间(秩3)', ...
           'F: STFT+Hough脊线+TF掩模', ...
           'D: STFT+Viterbi+TF掩模(本文)', ...
           'D-noIMU: Viterbi无IMU先验'};

for s = 1:n_schemes
    subplot(3,3,s);
    imagesc(range_axis, v_tf, 20*log10(abs(rd_list_viz{s})+eps)); axis xy; colorbar;
    clim(clim_rd); hold on;
    plot(range_axis(MID.r_ref), v_tf(MID.d_ref), 'r^','MarkerSize',10,'LineWidth',2);
    xlabel('距离 (m)'); ylabel('速度 (m/s)'); title(ttl_viz{s}); grid on;
end
subplot(3,3,8);
imagesc(range_axis, v_tf, 20*log10(abs(MID.rd_pre)+eps)); axis xy; colorbar;
clim(clim_rd); hold on;
plot(range_axis(MID.r_ref), v_tf(MID.d_ref), 'r^','MarkerSize',10,'LineWidth',2);
xlabel('距离 (m)'); ylabel('速度 (m/s)'); title('输入RD（未滤波）'); grid on;

subplot(3,3,9); axis off;
text(0.1,0.9,sprintf(['消融实验 — 第%d帧\n', ...
    '目标 R=%.0fm, v=%.1fm/s\n\n', ...
    'A-prime SCR=%.1f dB\n', ...
    'B-prime SCR=%.1f dB\n', ...
    'C:  SCR=%.1f dB\n', ...
    'E:  SCR=%.1f dB\n', ...
    'F:  SCR=%.1f dB\n', ...
    'D:  SCR=%.1f dB\n', ...
    'D-noIMU: SCR=%.1f dB'], ...
    f_mid, MID.tgt_R, MID.v_tgt, ...
    schemes.Ap.SCR(f_mid), schemes.Bp.SCR(f_mid), ...
    schemes.C.SCR(f_mid), schemes.E.SCR(f_mid), ...
    schemes.F.SCR(f_mid), schemes.D.SCR(f_mid), schemes.DnoIMU.SCR(f_mid)), ...
    'FontSize',10,'VerticalAlignment','top');
sgtitle(sprintf('七方案 RD 图消融对比（第%d帧, t=%.2fs）', f_mid, MID.t_fc), ...
    'FontSize',13,'FontWeight','bold');

% ── 图2：SCR 时序对比（7 方案）─────────────────────────────────────────────
figure('Name','SCR时序对比','Position',[40,80,1200,450]);
colors = lines(n_schemes);
for si = 1:n_schemes
    plot(t_frames, schemes.(scheme_names{si}).SCR, 'Color', colors(si,:), ...
        'LineWidth', 1.2+0.3*(si>=6), 'DisplayName', scheme_names{si}); hold on;
end
xlabel('时间 (s)'); ylabel('SCR (dB)'); title('各帧信杂比时序 — 7 方案消融对比');
legend('Location','best','FontSize',8); grid on;

% ── 图3：TERI 时序 ─────────────────────────────────────────────────────────
figure('Name','TERI时序','Position',[40,80,1200,450]);
for si = 3:n_schemes  % 仅自适应方案有 TERI 对比意义
    if si==1||si==2; continue; end  % A'/B' 硬零陷导致 TERI≈0
    plot(t_frames, schemes.(scheme_names{si}).TERI, 'Color', colors(si,:), ...
        'LineWidth', 1.5, 'DisplayName', scheme_names{si}); hold on;
end
yline(1.0,'k--','TERI=1.0','LineWidth',1.0);
xlabel('时间 (s)'); ylabel('目标能量保持率 TERI');
title('目标能量保持率 — 自适应方案对比');
legend('Location','best'); grid on;

% ── 图4：ROC 曲线（均值）───────────────────────────────────────────────────
figure('Name','ROC曲线','Position',[40,80,1000,450]);
subplot(1,2,1);  % Pd vs Pfa (全帧平均)
for si = 1:n_schemes
    sn = scheme_names{si};
    pd_mean = mean(roc_results.(sn).Pd_vs_Pfa, 1);
    plot(PFA_LIST, pd_mean, 'Color', colors(si,:), 'LineWidth', 1.2+0.3*(si>=6), ...
        'DisplayName', sn); hold on;
end
set(gca, 'XScale', 'log'); xlabel('Pfa (CFAR理论值)'); ylabel('检测概率 Pd');
title('ROC曲线 — Pd vs Pfa（全帧均值）'); legend('Location','best','FontSize',8);
grid on; xlim([min(PFA_LIST), max(PFA_LIST)]); ylim([0,1]);

subplot(1,2,2);  % 实际FAR vs Pd
for si = 1:n_schemes
    sn = scheme_names{si};
    far_mean = mean(roc_results.(sn).FAR_vs_Pfa, 1);
    pd_mean = mean(roc_results.(sn).Pd_vs_Pfa, 1);
    plot(far_mean, pd_mean, 'Color', colors(si,:), 'LineWidth', 1.2+0.3*(si>=6), ...
        'DisplayName', sn); hold on;
end
set(gca, 'XScale', 'log'); xlabel('实际虚警率 FAR'); ylabel('检测概率 Pd');
title('ROC曲线 — Pd vs 实际FAR（全帧均值）'); legend('Location','best','FontSize',8);
grid on; xlim([1e-5, 1]); ylim([0,1]);
sgtitle('消融实验 ROC 曲线评估（逐方案扫 CFAR Pfa）','FontSize',12,'FontWeight','bold');

% ── 图5：消融柱状图 ─────────────────────────────────────────────────────────
figure('Name','消融柱状图','Position',[40,80,1200,450]);
subplot(1,3,1);
bar(1:n_schemes, scr_avg, 0.6, 'FaceColor', [0.3 0.6 0.9]); hold on;
for si = 1:n_schemes
    text(si, scr_avg(si)+0.5, sprintf('%.1f', scr_avg(si)), ...
        'Hor','center','Vert','bottom','FontSize',8);
end
set(gca,'XTick',1:n_schemes,'XTickLabel',scheme_names,'XTickLabelRotation',30);
ylabel('平均SCR (dB)'); title('SCR对比'); grid on;
% 标注消融链箭头
annotation('textarrow',[0.25,0.38],[0.6,0.6],'String','RD域链','FontSize',8);
annotation('textarrow',[0.55,0.78],[0.6,0.6],'String','时频域链','FontSize',8);

subplot(1,3,2);
if nv>0
    bar(3:n_schemes, teri_avg(3:n_schemes), 0.6, 'FaceColor', [0.9 0.5 0.2]);
    set(gca,'XTick',3:n_schemes,'XTickLabel',scheme_names(3:n_schemes),'XTickLabelRotation',30);
    ylabel('TERI'); title('目标能量保持率'); grid on;
end

subplot(1,3,3);
bar(1:n_schemes, far_avg, 0.6, 'FaceColor', [0.7 0.2 0.2]); hold on;
set(gca,'XTick',1:n_schemes,'XTickLabel',scheme_names,'XTickLabelRotation',30);
ylabel('FAR'); title('虚警率对比'); grid on;
sgtitle('七方案消融实验汇总','FontSize',12,'FontWeight','bold');

fprintf('[完成] 所有图表已生成。\n');
fprintf('============================================================\n');
fprintf('  消融实验结论：\n');
fprintf('  (1) 时频面 vs RD面：  C vs D-noIMU： Δ=%.2f dB\n', scr_avg(7)-scr_avg(3));
fprintf('  (2) Viterbi vs Hough： F vs D-noIMU： Δ=%.2f dB\n', scr_avg(7)-scr_avg(5));
fprintf('  (3) IMU先验贡献：     D-noIMU vs D： Δ=%.2f dB\n', scr_avg(6)-scr_avg(7));
fprintf('============================================================\n');
end

%% =========================================================================
%  ====================  内 置 函 数 区  ===================================
%% =========================================================================

% ─────────────────────────────────────────────────────────────────────────
function rd_out = apply_conservative_notch(rd_in, v_axis, mu_c, B_core)
%APPLY_CONSERVATIVE_NOTCH  保守核心陷波（两次滤波第一阶段）
H_min_p1 = 0.01;
v_col    = v_axis(:);
H_vec    = ones(numel(v_col), 1);
core_idx = abs(v_col - mu_c) <= B_core;
H_vec(core_idx) = H_min_p1;
H        = repmat(H_vec, 1, size(rd_in, 2));
rd_out   = rd_in .* H;
end

% ─────────────────────────────────────────────────────────────────────────
function [rd_out, diag] = tf_ridge_filter(frm, win_doppler, tfc, dv_att_pulse, sigma_pred_v, tgt_range_bin)
%  Viterbi 脊线追踪 + 时频掩模（方案 D，同 v3）
%  完整实现见下方注释，与 v3 完全一致
[N_slow, N_fast] = size(frm);
L=tfc.L; H=tfc.H; Nfft=tfc.Nfft; pad=tfc.pad; w=tfc.w(:); v_tf=tfc.v_tf(:);
dv_bin = v_tf(2)-v_tf(1);

xp = [zeros(pad,N_fast); frm; zeros(pad,N_fast)];
Np = size(xp,1);
starts = 1:H:(Np-L+1);
Nf = numel(starts);
X_tf = zeros(Nfft, Nf, N_fast);
Ptf  = zeros(Nfft, Nf);
Ptf_tgt = zeros(Nfft, Nf);
if nargin>=6 && ~isempty(tgt_range_bin)
    tr = max(1,tgt_range_bin-1):min(N_fast,tgt_range_bin+1);
else
    tr = [];
end
for i = 1:Nf
    s = starts(i);
    seg = xp(s:s+L-1,:) .* w;
    F   = fftshift(fft(seg, Nfft, 1), 1);
    X_tf(:,i,:) = F;
    Ptf(:,i) = sum(abs(F(:,tfc.clut_bins)).^2, 2);
    if ~isempty(tr)
        Ptf_tgt(:,i) = sum(abs(F(:,tr)).^2, 2);
    end
end

% 脊线预测：数据锚点 + IMU 蜿蜒先验
Ptf_sum = sum(Ptf,2);
[~, k_anchor] = max(Ptf_sum);
delta_v = zeros(Nf,1);
for n = 1:Nf
    c0 = starts(n) - pad + L/2;
    a = min(max(round(c0-L/2),1),N_slow);
    b = min(max(round(c0+L/2),1),N_slow);
    if b<a; tmp=a; a=b; b=tmp; end
    delta_v(n) = mean(dv_att_pulse(a:b));
end
delta_v = delta_v - mean(delta_v);
k_pred  = round(k_anchor + delta_v/dv_bin);
k_pred  = min(max(k_pred,1),Nfft);

% Viterbi 脊线追踪
k_ridge = viterbi_ridge(Ptf, k_pred, tfc.dk_search, tfc.sigmaT, tfc.lamT);

% 逐帧 3dB 谱宽 + 融合 → 平顶/滚降半宽
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

% 掩模 G(k,n)
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

% 距离自适应 CNR 权重 w_r
in_band_kv  = abs((1:Nfft).' - round(mean(k_ridge))) <= round(mean(B_flat_bin));
out_band_kv = ~in_band_kv;
n_in = sum(in_band_kv); n_out = max(sum(out_band_kv),1);
P_in_r  = zeros(1, N_fast);
P_out_r = zeros(1, N_fast);
for i = 1:Nf
    F_i = squeeze(X_tf(:,i,:));
    P_in_r  = P_in_r  + sum(abs(F_i(in_band_kv,:)).^2, 1);
    P_out_r = P_out_r + sum(abs(F_i(out_band_kv,:)).^2, 1);
end
P_out_norm   = P_out_r * (n_in/n_out);
CNR_excess   = max(0, P_in_r - P_out_norm) ./ (P_out_norm + eps);
w_r          = CNR_excess ./ (1 + CNR_excess);
w_r          = movmean(w_r, 7);
w_r          = min(max(w_r, 0), 1);

% WOLA-ISTFT 重建
yp = zeros(Np, N_fast); ws = zeros(Np,1);
for i = 1:Nf
    s = starts(i);
    F_i  = squeeze(X_tf(:,i,:));
    G_r  = 1 + w_r .* (G(:,i) - 1);
    col  = F_i .* G_r;
    seg  = ifft(ifftshift(col, 1), Nfft, 1);
    seg  = seg(1:L,:);
    yp(s:s+L-1,:) = yp(s:s+L-1,:) + seg .* w;
    ws(s:s+L-1)   = ws(s:s+L-1) + w.^2;
end
ws(ws<1e-12) = 1e-12;
frm_filt = yp(pad+1:pad+N_slow,:) ./ ws(pad+1:pad+N_slow);

% 多普勒 FFT → RD
rd_out = fftshift(fft(frm_filt .* win_doppler, N_slow, 1), 1);

diag.Ptf=Ptf; diag.Ptf_tgt=Ptf_tgt; diag.Ptf_masked=Ptf.*G.^2;
diag.k_pred=k_pred; diag.k_ridge=k_ridge; diag.G=G; diag.Nf=Nf;
diag.v_ridge = v_tf(k_ridge); diag.B_flat_v = B_flat_bin*dv_bin; diag.B_roll_v = B_roll_bin*dv_bin;
diag.w_r = w_r; diag.CNR_excess = CNR_excess;
diag.B_flat_used_ms = mean(B_flat_bin) * dv_bin;
diag.B_flat_max_ms  = max(B_flat_bin)  * dv_bin;
end

% ─────────────────────────────────────────────────────────────────────────
function [rd_out, diag] = tf_ridge_filter_noimu(frm, win_doppler, tfc, sigma_pred_v, tgt_range_bin)
%  与 tf_ridge_filter 相同，但 Viterbi 中去除 IMU 先验（方案 D-noIMU）
%  关键区别：k_pred 固定为 k_anchor（纯数据驱动），不引入逐帧 delta_v
[N_slow, N_fast] = size(frm);
L=tfc.L; H=tfc.H; Nfft=tfc.Nfft; pad=tfc.pad; w=tfc.w(:); v_tf=tfc.v_tf(:);
dv_bin = v_tf(2)-v_tf(1);

xp = [zeros(pad,N_fast); frm; zeros(pad,N_fast)];
Np = size(xp,1);
starts = 1:H:(Np-L+1);
Nf = numel(starts);
X_tf = zeros(Nfft, Nf, N_fast);
Ptf  = zeros(Nfft, Nf);
Ptf_tgt = zeros(Nfft, Nf);
if nargin>=6 && ~isempty(tgt_range_bin)
    tr = max(1,tgt_range_bin-1):min(N_fast,tgt_range_bin+1);
else
    tr = [];
end
for i = 1:Nf
    s = starts(i);
    seg = xp(s:s+L-1,:) .* w;
    F   = fftshift(fft(seg, Nfft, 1), 1);
    X_tf(:,i,:) = F;
    Ptf(:,i) = sum(abs(F(:,tfc.clut_bins)).^2, 2);
    if ~isempty(tr)
        Ptf_tgt(:,i) = sum(abs(F(:,tr)).^2, 2);
    end
end

% ── 关键区别：无 IMU 先验，所有帧 k_pred = k_anchor ──
Ptf_sum = sum(Ptf,2);
[~, k_anchor] = max(Ptf_sum);
k_pred = round(k_anchor) * ones(Nf, 1);  % 纯数据驱动，无蜿蜒先验
k_pred = min(max(k_pred,1),Nfft);

% Viterbi（使用常值先验）
k_ridge = viterbi_ridge(Ptf, k_pred, tfc.dk_search, tfc.sigmaT, tfc.lamT);

% 以下与 tf_ridge_filter 完全相同
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

in_band_kv  = abs((1:Nfft).' - round(mean(k_ridge))) <= round(mean(B_flat_bin));
out_band_kv = ~in_band_kv;
n_in = sum(in_band_kv); n_out = max(sum(out_band_kv),1);
P_in_r  = zeros(1, N_fast);
P_out_r = zeros(1, N_fast);
for i = 1:Nf
    F_i = squeeze(X_tf(:,i,:));
    P_in_r  = P_in_r  + sum(abs(F_i(in_band_kv,:)).^2, 1);
    P_out_r = P_out_r + sum(abs(F_i(out_band_kv,:)).^2, 1);
end
P_out_norm   = P_out_r * (n_in/n_out);
CNR_excess   = max(0, P_in_r - P_out_norm) ./ (P_out_norm + eps);
w_r          = CNR_excess ./ (1 + CNR_excess);
w_r          = movmean(w_r, 7);
w_r          = min(max(w_r, 0), 1);

yp = zeros(Np, N_fast); ws = zeros(Np,1);
for i = 1:Nf
    s = starts(i);
    F_i  = squeeze(X_tf(:,i,:));
    G_r  = 1 + w_r .* (G(:,i) - 1);
    col  = F_i .* G_r;
    seg  = ifft(ifftshift(col, 1), Nfft, 1);
    seg  = seg(1:L,:);
    yp(s:s+L-1,:) = yp(s:s+L-1,:) + seg .* w;
    ws(s:s+L-1)   = ws(s:s+L-1) + w.^2;
end
ws(ws<1e-12) = 1e-12;
frm_filt = yp(pad+1:pad+N_slow,:) ./ ws(pad+1:pad+N_slow);

rd_out = fftshift(fft(frm_filt .* win_doppler, N_slow, 1), 1);

diag.Ptf=Ptf; diag.Ptf_tgt=Ptf_tgt; diag.Ptf_masked=Ptf.*G.^2;
diag.k_pred=k_pred; diag.k_ridge=k_ridge; diag.G=G; diag.Nf=Nf;
diag.v_ridge = v_tf(k_ridge); diag.B_flat_v = B_flat_bin*dv_bin; diag.B_roll_v = B_roll_bin*dv_bin;
diag.w_r = w_r; diag.CNR_excess = CNR_excess;
diag.B_flat_used_ms = mean(B_flat_bin) * dv_bin;
diag.B_flat_max_ms  = max(B_flat_bin)  * dv_bin;
end

% ─────────────────────────────────────────────────────────────────────────
function rd_out = pca_clutter_suppression(frm, win_doppler, n_components, N_slow,N_fast)
%PCA_CLUTTER_SUPPRESSION  PCA 子空间投影杂波抑制（方案 E）
%  输入 frm: N_slow × N_fast 慢时间-距离矩阵（MOCO 后）
%  原理：对慢时间协方差矩阵做特征分解，前 n_components 个主成分对应杂波子空间，
%        正交投影抑制杂波。
%
%  实现步骤：
%    (1) 对每个距离门，提取慢时间序列
%    (2) 构建数据矩阵 X [N_slow × N_fast]
%    (3) 估计协方差矩阵 R = (1/N_fast) * X * X^H
%    (4) 特征分解，取前 n_components 个特征向量为杂波子空间 U_c
%    (5) 投影矩阵 P_c = U_c * U_c^H
%    (6) 输出 Y = X - P_c * X（抑制杂波）
%
%  注：PCA 不依赖 IMU 先验，是公平的非参数对比方法

X = frm.';                        % [N_fast × N_slow]
X = X - mean(X, 1);               % 去直流
R = (X' * X) / N_fast;            % [N_slow × N_slow] 协方差矩阵

[U, ~, ~] = svd(R);
n_comp = min(n_components, N_slow-1);
U_c = U(:, 1:n_comp);             % 杂波子空间
P_c = U_c * U_c';                 % 投影矩阵

Y = X - X * P_c;                  % 正交投影抑制
frm_out = Y.';                    % 转回 [N_slow × N_fast]

% 加窗 → FFT → RD
rd_out = fftshift(fft(frm_out .* win_doppler, N_slow, 1), 1);
end

% ─────────────────────────────────────────────────────────────────────────
function [rd_out, diag] = hough_ridge_filter(frm, win_doppler, tfc, tgt_range_bin, houghOpts)
%HOUGH_RIDGE_FILTER  Hough 变换脊线追踪 + 时频掩模（方案 F）
%  与 tf_ridge_filter 相同框架，但用 Hough 变换替代 Viterbi 做脊线检测
%
%  原理：
%    (1) 对时频谱 Ptf 做二值化（阈值为 energy_th × 峰值）
%    (2) 对二值化图像做 Hough 变换（MATLAB hough()）
%    (3) 取峰值对应的直线参数 (θ, ρ) → 映射为帧索引 k_ridge
%    (4) 以此脊线构造时频掩模 G(k,n)
%    (5) WOLA-ISTFT 重建
%
%  与 Viterbi 的关键区别：
%      Hough 只能检测全局直线（或通过窗内分段近似曲线），
%      Viterbi 可追踪任意蜿蜒路径且受 IMU 先验约束。
%      通过 F vs D-noIMU 对比，可分离"全局 Hough"和"局部 Viterbi"的差异。

[N_slow, N_fast] = size(frm);
L=tfc.L; H=tfc.H; Nfft=tfc.Nfft; pad=tfc.pad; w=tfc.w(:); v_tf=tfc.v_tf(:);
dv_bin = v_tf(2)-v_tf(1);

% STFT（同 tf_ridge_filter）
xp = [zeros(pad,N_fast); frm; zeros(pad,N_fast)];
Np = size(xp,1);
starts = 1:H:(Np-L+1);
Nf = numel(starts);
X_tf = zeros(Nfft, Nf, N_fast);
Ptf  = zeros(Nfft, Nf);
Ptf_tgt = zeros(Nfft, Nf);
if nargin>=4 && ~isempty(tgt_range_bin)
    tr = max(1,tgt_range_bin-1):min(N_fast,tgt_range_bin+1);
else
    tr = [];
end
for i = 1:Nf
    s = starts(i);
    seg = xp(s:s+L-1,:) .* w;
    F   = fftshift(fft(seg, Nfft, 1), 1);
    X_tf(:,i,:) = F;
    Ptf(:,i) = sum(abs(F(:,tfc.clut_bins)).^2, 2);
    if ~isempty(tr)
        Ptf_tgt(:,i) = sum(abs(F(:,tr)).^2, 2);
    end
end

% ── Hough 变换脊线检测 ──
% 将 Ptf 归一化后二值化
Ptf_norm = Ptf / (max(Ptf(:)) + eps);
bw_img = Ptf_norm > houghOpts.energy_th;

% Hough 变换（theta范围: [-90, 89]，MATLAB hough不接受90）
theta_range = -90 + (0:houghOpts.n_theta-1) * 180/houghOpts.n_theta;
[H_mat, theta_h, rho_h] = hough(bw_img, ...
    'Theta', theta_range, ...
    'RhoResolution', (size(bw_img,1)+size(bw_img,2))/houghOpts.n_rho);

peaks = houghpeaks(H_mat, houghOpts.npeaks, 'Threshold', 0.5*max(H_mat(:)));
if isempty(peaks)
    % Hough 未检测到脊线，回退到能量峰值
    [~, k_ridge_raw] = max(Ptf, [], 1);
    k_ridge = round(mean(k_ridge_raw)) * ones(Nf, 1);
else
    lines_h = houghlines(bw_img, theta_h, rho_h, peaks, 'FillGap', 5, 'MinLength', 3);
    if isempty(lines_h)
        [~, k_ridge_raw] = max(Ptf, [], 1);
        k_ridge = round(mean(k_ridge_raw)) * ones(Nf, 1);
    else
        % 取最长线段对应的直线参数
        len_max = 0; best_line = [];
        for li = 1:length(lines_h)
            % 从端点计算线段长度（houghlines 不直接提供长字段）
            dx = lines_h(li).point2(1) - lines_h(li).point1(1);
            dy = lines_h(li).point2(2) - lines_h(li).point1(2);
            seg_len = sqrt(dx^2 + dy^2);
            if seg_len > len_max
                len_max = seg_len;
                best_line = lines_h(li);
            end
        end
        % 从 Hough 直线参数 (θ, ρ) 计算每帧的 k_ridge
        theta_r = deg2rad(best_line.theta);
        k_ridge = zeros(Nf,1);
        for n = 1:Nf
            rho = best_line.rho;
            % 直线方程: x*cos(θ) + y*sin(θ) = ρ
            % x = n (帧索引), y = k (频率 bin)
            % k = (ρ - n*cos(θ)) / sin(θ)
            if abs(sin(theta_r)) < 1e-6
                k_ridge(n) = round(rho / cos(theta_r));
            else
                k_ridge(n) = round((rho - n*cos(theta_r)) / sin(theta_r));
            end
        end
        k_ridge = min(max(k_ridge, 1), Nfft);
    end
end

% ── 逐帧谱宽估计（同 tf_ridge_filter）──
% sigma_pred_bin = max(1, 2.0);  % Hough 无双 IMU 先验，用保守值
B_flat_bin = zeros(Nf,1);
B_roll_bin = zeros(Nf,1);
for n = 1:Nf
    col = Ptf(:,n);
    pk = col(k_ridge(n))+eps;
    above = col >= 0.1*pk;
    bw = local_3db_width(above, k_ridge(n));
    sig_data_bin = max(bw/1.5, 2.0);
    B_flat_bin(n) = min(tfc.kappa*sig_data_bin, tfc.Bflat_cap/dv_bin);
    B_roll_bin(n) = 0.5*sig_data_bin;
end

% ── 掩模 G(k,n) ──
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

% ── 距离自适应 CNR 权重 w_r ──
in_band_kv  = abs((1:Nfft).' - round(mean(k_ridge))) <= round(mean(B_flat_bin));
out_band_kv = ~in_band_kv;
n_in = sum(in_band_kv); n_out = max(sum(out_band_kv),1);
P_in_r  = zeros(1, N_fast);
P_out_r = zeros(1, N_fast);
for i = 1:Nf
    F_i = squeeze(X_tf(:,i,:));
    P_in_r  = P_in_r  + sum(abs(F_i(in_band_kv,:)).^2, 1);
    P_out_r = P_out_r + sum(abs(F_i(out_band_kv,:)).^2, 1);
end
P_out_norm   = P_out_r * (n_in/n_out);
CNR_excess   = max(0, P_in_r - P_out_norm) ./ (P_out_norm + eps);
w_r          = CNR_excess ./ (1 + CNR_excess);
w_r          = movmean(w_r, 7);
w_r          = min(max(w_r, 0), 1);

% ── WOLA-ISTFT 重建 ──
yp = zeros(Np, N_fast); ws = zeros(Np,1);
for i = 1:Nf
    s = starts(i);
    F_i  = squeeze(X_tf(:,i,:));
    G_r  = 1 + w_r .* (G(:,i) - 1);
    col  = F_i .* G_r;
    seg  = ifft(ifftshift(col, 1), Nfft, 1);
    seg  = seg(1:L,:);
    yp(s:s+L-1,:) = yp(s:s+L-1,:) + seg .* w;
    ws(s:s+L-1)   = ws(s:s+L-1) + w.^2;
end
ws(ws<1e-12) = 1e-12;
frm_filt = yp(pad+1:pad+N_slow,:) ./ ws(pad+1:pad+N_slow);

rd_out = fftshift(fft(frm_filt .* win_doppler, N_slow, 1), 1);

diag.Ptf=Ptf; diag.Ptf_tgt=Ptf_tgt; diag.Ptf_masked=Ptf.*G.^2;
diag.k_ridge=k_ridge; diag.G=G; diag.Nf=Nf;
diag.v_ridge = v_tf(k_ridge); diag.B_flat_v = B_flat_bin*dv_bin;
diag.w_r = w_r; diag.CNR_excess = CNR_excess;
diag.B_flat_used_ms = mean(B_flat_bin) * dv_bin;
end

% ─────────────────────────────────────────────────────────────────────────
function ridge = viterbi_ridge(Ptf, k_pred, dk, sigmaT, lamT)
%VITERBI_RIDGE  IMU 先验约束 + 能量驱动的时频脊线动态规划追踪（同 v3）
[Kf, Nf] = size(Ptf);
kk = (1:Kf).';
sigp = max(dk/2, 1);
E = log(Ptf + 1e-12) - ((kk - k_pred(:).').^2)/(2*sigp^2);
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
if ~above(kc); bw = 1; return; end
lo = kc;
while lo>1 && above(lo-1); lo=lo-1; end
hi = kc;
while hi<N && above(hi+1); hi=hi+1; end
bw = hi - lo + 1;
end

% ─────────────────────────────────────────────────────────────────────────
function H = adaptive_soft_notch(v_axis, mu_c, B_flat, B_roll, alpha, H_min, N_dop, N_rng)
%ADAPTIVE_SOFT_NOTCH  目标失真约束的自适应软陷波（平顶+高斯滚降）
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
%SIM_DYNAMIC_CLUTTER  物理散射单元模型的动态地面杂波（同 v3）
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
%  从 RD 图估计残余杂波中心与谱宽（同 v3）
if nargin<3 || isempty(mu_prior);  mu_prior = 0;   end
if nargin<4 || isempty(search_hw); search_hw = 4.0; end
v_axis = v_axis(:);
if nargin>=5 && ~isempty(clut_bins)
    cb  = clut_bins(clut_bins>=1 & clut_bins<=size(rd,2));
    if isempty(cb); cb = 1:size(rd,2); end
    pwr = mean(abs(rd(:, cb)).^2, 2);
else
    pwr = mean(abs(rd).^2, 2);
end
in_win = abs(v_axis - mu_prior) <= search_hw;
if nnz(in_win) < 5; mu_c = mu_prior; sigma_c = 1.0; return; end
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
%COMPUTE_CLUTTER_RIDGE_MODEL  基于 IMU 运动状态的每 CPI 脊线参数（同 v3）
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
    SIGMA_C_FLOOR = 0.8;
    sigma_c = max(sqrt(sigma_c_beam^2 + sigma_c_cpi^2), SIGMA_C_FLOOR);
    model.sigma_c_est(f) = sigma_c;
    model.clutter_pwr_ratio(f) = min(1.0, v_mag_plat/15);
end
end

% ─────────────────────────────────────────────────────────────────────────
function [d_found, r_found, d_refined, r_refined, found_ok] = find_target_peak_prior(rd, d_prior, r_prior, W_d, W_r, gamma_min_dB)
%FIND_TARGET_PEAK_PRIOR  先验引导的目标峰值搜索（同 v3）
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
function [scr_dB, P_tgt] = compute_scr(rd, tgt_dop, tgt_bin, N_dop, N_rng)
%COMPUTE_SCR  信杂比 (dB)（同 v3）
rb=3; db=3;
r_idx = max(1,tgt_bin-rb):min(N_rng,tgt_bin+rb);
d_idx = max(1,tgt_dop-db):min(N_dop,tgt_dop+db);
P_tgt = max(abs(rd(d_idx,r_idx)).^2, [], 'all');
P_tgt = max(P_tgt, mean(abs(rd(d_idx,r_idx)).^2,'all'));
dc = round(N_dop/2); cld = max(1,dc-25):min(N_dop,dc+25);
P_c = mean(abs(rd(cld,:)).^2,'all');
scr_dB = 10*log10(P_tgt/(P_c+eps));
end

% ─────────────────────────────────────────────────────────────────────────
function bw = compute_bw3dB(rd, clut_bins, v_axis, N_dop)
%COMPUTE_BW3DB  杂波多普勒 -3dB 带宽（m/s）（同 v3）
pwr = zeros(N_dop,1);
for b = clut_bins
    if b <= size(rd,2); pwr = pwr + abs(rd(:,b)).^2; end
end
pwr_n = pwr/(max(pwr)+eps);
bw = sum(pwr_n>=0.1)*(v_axis(2)-v_axis(1));
end

% ─────────────────────────────────────────────────────────────────────────
function det = my_go_cfar(pwr_map, Tr, Gr, Td, Gd, Pfa)
%MY_GO_CFAR  GO-CFAR 检测（同 v3）
[N_dop, N_rng] = size(pwr_map);
det = false(N_dop, N_rng);
N_ref = 2*Tr*Td; alpha_cfar = N_ref*(Pfa^(-1/N_ref)-1);
for di = Td+Gd+1 : N_dop-Td-Gd
    for ri = Tr+Gr+1 : N_rng-Tr-Gr
        d_lo_s=di-Td-Gd; d_lo_e=di-Gd-1; d_hi_s=di+Gd+1; d_hi_e=di+Gd+Td;
        r_lo_s=ri-Tr-Gr; r_lo_e=ri-Gr-1; r_hi_s=ri+Gr+1; r_hi_e=ri+Gr+Tr;
        if d_lo_s<1||d_hi_e>N_dop||r_lo_s<1||r_hi_e>N_rng; continue; end
        Z_d = sum(pwr_map(d_lo_s:d_lo_e,ri)) + sum(pwr_map(d_hi_s:d_hi_e,ri));
        Z_r = sum(pwr_map(di,r_lo_s:r_lo_e)) + sum(pwr_map(di,r_hi_s:r_hi_e));
        Z = max(max(Z_d,Z_r), eps);
        if pwr_map(di,ri) > alpha_cfar*Z/N_ref; det(di,ri)=true; end
    end
end
end

% ─────────────────────────────────────────────────────────────────────────
function [det_ref, totals_pfa, pd_pfa] = batch_go_cfar_roc(pwr_map, Tr, Gr, Td, Gd, Pfa_ref, Pfa_list, d_ref, r_ref)
%BATCH_GO_CFAR_ROC  单次遍历同时计算参考Pfa检测图和ROC多Pfa统计量
%  相比循环调用 my_go_cfar 21次，本函数只遍历数据一次，节约约 20× 开销。
%
%  输入：
%    pwr_map   — 功率图 [N_dop × N_rng]
%    Tr,Gr,Td,Gd — CFAR 窗参数
%    Pfa_ref   — 参考检测门限（输出完整检测图用于 pts/Pd 统计）
%    Pfa_list  — [1×N_pfa] ROC 扫描 Pfa 值
%    d_ref,r_ref — 目标单元索引（仅在该位置记录 Pd）
%  输出：
%    det_ref   — Pfa_ref 对应的完整检测图
%    totals_pfa — [1×N_pfa] 各 Pfa 总检测数
%    pd_pfa    — [1×N_pfa] 目标单元在各 Pfa 是否被检测到

N_pfa = length(Pfa_list);
N_ref = 2*Tr*Td;
alpha_ref = N_ref * (Pfa_ref^(-1/N_ref) - 1);
alpha_list = N_ref * (Pfa_list.^(-1/N_ref) - 1);

[N_dop, N_rng] = size(pwr_map);
det_ref = false(N_dop, N_rng);
totals_pfa = zeros(1, N_pfa);
pd_pfa = false(1, N_pfa);

for di = Td+Gd+1 : N_dop-Td-Gd
    for ri = Tr+Gr+1 : N_rng-Tr-Gr
        d_lo_s=di-Td-Gd; d_lo_e=di-Gd-1; d_hi_s=di+Gd+1; d_hi_e=di+Gd+Td;
        r_lo_s=ri-Tr-Gr; r_lo_e=ri-Gr-1; r_hi_s=ri+Gr+1; r_hi_e=ri+Gr+Tr;
        if d_lo_s<1||d_hi_e>N_dop||r_lo_s<1||r_hi_e>N_rng; continue; end
        Z_d = sum(pwr_map(d_lo_s:d_lo_e,ri)) + sum(pwr_map(d_hi_s:d_hi_e,ri));
        Z_r = sum(pwr_map(di,r_lo_s:r_lo_e)) + sum(pwr_map(di,r_hi_s:r_hi_e));
        Z = max(max(Z_d,Z_r), eps) / N_ref;

        val = pwr_map(di,ri);

        % 参考 Pfa：记录完整检测图
        if val > alpha_ref * Z
            det_ref(di,ri) = true;
        end

        % ROC 各 Pfa：只累计统计量，不存全图
        is_tgt = (di == d_ref) && (ri == r_ref);
        for pi = 1:N_pfa
            if val > alpha_list(pi) * Z
                totals_pfa(pi) = totals_pfa(pi) + 1;
                if is_tgt; pd_pfa(pi) = true; end
            end
        end
    end
end
end

% ─────────────────────────────────────────────────────────────────────────
function [cand_mask, cand_list] = detect_candidate_targets(rd_power, thresh_db)
%DETECT_CANDIDATE_TARGETS  局部峰值检测 + 非极大值抑制（同 v3）
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
%ISOLATED_PEAK_TEST  孤立点目标判定（同 v3）
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
