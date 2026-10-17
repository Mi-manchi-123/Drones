function results = run_single_experiment(SNR_db, mc_seed)
%RUN_SINGLE_EXPERIMENT  单次蒙特卡洛实验
%   运行四方案(A/B/C/D)全部处理流程，返回性能指标
%   输入:
%     SNR_db  — 信噪比 (dB)
%     mc_seed — 随机种子（控制噪声实现，杂波种子固定在 42）
%   输出:
%     results — 结构体包含 .SCR_mean, .Pd, .FAR, .TERI, .AUC, .Pd_roc, .FAR_roc 等

%% =========================================================================
% PART 0：公共雷达参数
% =========================================================================
Ob_time  = 2.0;
fc       = 13.5e9;
B        = 20e6;
Tc       = 5e-6;
PRI      = 200e-6;
fs       = 60e6;
c        = 3e8;
lambda   = c / fc;
K        = B / Tc;

N_fast   = 1024;
N_slow   = 256;
overlap  = 100;
Num_chirp = floor(Ob_time / PRI);

t_fast      = (0:N_fast-1) / fs;
range_axis  = t_fast * c / 2;
range_res   = c / (fs * 2);

step       = N_slow - overlap;
num_frames = floor((Num_chirp - N_slow) / step) + 1;
v_tf       = ((0:N_slow-1) - floor(N_slow/2)) * (1/PRI) / N_slow * lambda/2;
v_tf       = v_tf(:).';

win_range   = hanning(N_fast).';
win_doppler = hanning(N_slow);

%% =========================================================================
% PART 1：平台运动状态生成
% =========================================================================
t_sim    = (0:Num_chirp-1).' * PRI;
v_plat   = 10;
heading  = 30;
P_x = v_plat * cos(deg2rad(heading)) * t_sim;
P_y = v_plat * sin(deg2rad(heading)) * t_sim;
P_z = zeros(Num_chirp, 1);
A_x = 3.0 * sin(2*pi*1.5*t_sim) + 0.2*sin(2*pi*25*t_sim + 0.7);
A_y = 1.5 * sin(2*pi*1.2*t_sim + 1.0) + 0.1*sin(2*pi*30*t_sim);
A_z = heading + 0.5*sin(2*pi*0.5*t_sim);
eskf_data.t   = t_sim;
eskf_data.P_x = P_x; eskf_data.P_y = P_y; eskf_data.P_z = P_z;
eskf_data.A_x = A_x; eskf_data.A_y = A_y; eskf_data.A_z = A_z;

%% =========================================================================
% PART 2：目标参数
% =========================================================================
tgt(1).coor0 = [300, 200, 100];
tgt(1).vt    = [3.2, 2.4, 0];
tgt(1).A     = [0, 0, 0];
tgt(1).type  = 'const';
N_tgt        = 1;
radar_offset = [0, 0, 5];

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
IN_coor0 = tgt(1).coor0;
tgt_R0   = norm(tgt(1).coor0 - ant_pos_all(1,:));

s_rx_raw    = zeros(Num_chirp, N_fast);
target_pos  = zeros(Num_chirp, 3);
dist_target_arr = zeros(Num_chirp, 1);

for i = 1:Num_chirp
    t_curr       = (i-1) * PRI;
    ant_pos      = ant_pos_all(i,:);
    tgt_pos_i    = tgt(1).coor0 + t_curr * tgt(1).vt + 0.5 * tgt(1).A * t_curr^2;
    target_pos(i,:) = tgt_pos_i;
    dist_t       = norm(tgt_pos_i - ant_pos);
    dist_target_arr(i) = dist_t;

    range_amp = (tgt_R0 / dist_t)^2;
    tau       = 2 * dist_t / c;
    s_rx_raw(i,:) = s_rx_raw(i,:) + ...
        range_amp * rectpuls(t_fast - Tc/2 - tau, Tc) .* ...
        exp(1j*pi*K*(t_fast - tau).^2) .* exp(-1j*2*pi*fc*tau);
end

%% =========================================================================
% PART 5：生成/加载杂波（缓存 .mat 避免重复生成）
% =========================================================================
CLUTTER_CACHE = 's_clutter_full.mat';
if exist(CLUTTER_CACHE, 'file')
    load(CLUTTER_CACHE, 's_clutter_full');
else
    rng(42);
    s_clutter_full = sim_dynamic_clutter(Num_chirp, N_fast, fs, fc, PRI, K, Tc, ...
        t_fast, ant_pos_all, ant_vel_all, lambda, ...
        tgt(1).coor0, -10);
    save(CLUTTER_CACHE, 's_clutter_full');
end

%% =========================================================================
% PART 6：运动补偿（MOCO）+ 脉冲压缩（使用 mc_seed 控制噪声）
% =========================================================================
dist1 = vecnorm(ant_pos_all  - target_pos, 2, 2);
dist2 = vecnorm(ant_ideal_all - target_pos, 2, 2);
delta_swing_all     = dist1 - dist2;
swing_phase_2d_corr = exp(1j * 4*pi * delta_swing_all / lambda);

t_ref = (0:N_fast-1) / fs;
s_t   = rectpuls(t_ref - Tc/2, Tc) .* exp(1j*pi*K*t_ref.^2);
X_tx  = fft(s_t, N_fast);

% ---- MC 关键：噪声种子由 mc_seed 控制 ----
rng(mc_seed);
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

%% =========================================================================
% PART 7：运动状态驱动的杂波脊线参数
% =========================================================================
ant0_approx = pos_all(1,:) + radar_offset;
beam_dir    = atan2(IN_coor0(2)-ant0_approx(2), IN_coor0(1)-ant0_approx(1));
theta_bw    = deg2rad(30);
clutter_model = compute_clutter_ridge_model(vel_all, ant_vel_all, ang_rate_all, ...
    pos_all, ang_all, lambda, PRI, N_slow, step, num_frames, Num_chirp, radar_offset, ...
    beam_dir, theta_bw, IN_coor0, 0);

%% =========================================================================
% PART 8：处理配置
% =========================================================================
BETA_FUSION = 0.6;
KAPPA_BW    = 2.5;
ALPHA_EMA      = 0.7;
sigma_c_ema    = 1.5;
sigma_c_total_arr = zeros(num_frames,1);

tfc.L         = 64;
tfc.H         = 16;
tfc.Nfft      = N_slow;
tfc.pad       = tfc.L/2;
tfc.w         = hann(tfc.L);
tfc.v_tf      = v_tf;
tfc.dk_search = 30;
tfc.sigmaT    = 3;
tfc.lamT      = 0.4;
tfc.kappa     = KAPPA_BW;
tfc.Hmin      = 0.01;
tfc.alpha     = 1.0;
tfc.Bflat_cap = 3;
BETA_CONSERVATIVE = 0.8;
DELTA_MARGIN = 0.1;

%% =========================================================================
% PART 9：主处理循环
% =========================================================================
SCR_A = zeros(num_frames,1); SCR_B = zeros(num_frames,1);
SCR_C = zeros(num_frames,1); SCR_D = zeros(num_frames,1);
TERI_C = zeros(num_frames,1); TERI_D = zeros(num_frames,1);
TERI_A = zeros(num_frames,1); TERI_B = zeros(num_frames,1);
tgt_separable = false(num_frames,1);
pts_A = zeros(num_frames,1); pts_B = zeros(num_frames,1);
pts_C = zeros(num_frames,1); pts_D = zeros(num_frames,1);
BW3_C = zeros(num_frames,1); BW3_D = zeros(num_frames,1);
mu_c_all = zeros(num_frames,1);
sigma_c_all = zeros(num_frames,1);
peak_found = false(num_frames,1);

W_d_search = 4; W_r_search = 5; gamma_min_dB = 4.0;
Tr_go=8; Gr_go=2; Td_go=4; Gd_go=2; Pfa_work=1e-4;
WIDE_NOTCH_BW = 3.0;

% ROC sweeep init
PFA_SWEEP = logspace(-6, -1.5, 25);
nPfa = numel(PFA_SWEEP);
N_ref_cfar = 2*max(Td_go, Tr_go);
d_cfar_lo = Td_go+Gd_go+1; d_cfar_hi = N_slow-Td_go-Gd_go;
r_cfar_lo = Tr_go+Gr_go+1; r_cfar_hi = N_fast-Tr_go-Gr_go;
n_valid_cfar = (d_cfar_hi-d_cfar_lo+1)*(r_cfar_hi-r_cfar_lo+1);
valid_cfar_mask = false(N_slow, N_fast);
valid_cfar_mask(d_cfar_lo:d_cfar_hi, r_cfar_lo:r_cfar_hi) = true;
pd_count = zeros(nPfa, 4);
fa_roc_cnt = zeros(nPfa, 4);
pd_sep_cnt = zeros(nPfa, 4);
pd_clut_cnt = zeros(nPfa, 4);
n_sep_frames = 0;
n_clut_frames = 0;

for f_idx = 1:num_frames
    s_idx = (f_idx-1)*step + 1;
    e_idx = s_idx + N_slow - 1;
    if e_idx > Num_chirp; break; end
    frame_idx = s_idx:e_idx;

    m_fc = min(max(round((s_idx+e_idx)/2),1), Num_chirp);
    t_fc = (m_fc-1)*PRI;
    R_fc = build_R_zyx_local(deg2rad(ang_all(m_fc,1)), deg2rad(ang_all(m_fc,2)), deg2rad(ang_all(m_fc,3)));
    ant_pos_fc = pos_all(m_fc,:) + (R_fc*radar_offset.').';
    slow_idx = (0:N_slow-1).';

    e_los_bs   = R_fc;
    v_rad_plat = mean(vel_all(frame_idx,:) * e_los_bs, 1);
    fd_shift   = 2 * norm(v_rad_plat) / lambda;
    shift_plat = exp(-1j * 2*pi * fd_shift * slow_idx * PRI);

    frm_raw  = s_range_raw(frame_idx,:);
    frm_moco = s_range_moco(frame_idx,:);
    frm_raw_comp  = frm_raw  .* shift_plat;
    frm_moco_comp = frm_moco .* shift_plat;

    rd_pre = fftshift(fft(frm_moco_comp .* win_doppler, N_slow, 1), 1);

    mu_c_att_prior = clutter_model.mu_c_att(f_idx);
    clut_lo_est = 50;
    clut_hi_est = min(N_fast, round(0.6*N_fast));
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

    v_mean_f     = mean(vel_all(frame_idx,:),1);
    e_bore       = v_mean_f/(norm(v_mean_f)).';
    v_arm_f      = ant_vel_all(frame_idx,:) - vel_all(frame_idx,:);
    dv_att_pulse = -(v_arm_f .* e_bore);

    sigma_c_platform = clutter_model.sigma_c_est(f_idx);
    sigma_c_meas     = max(sigma_c_platform, sigma_c_data);
    sigma_c_ema      = ALPHA_EMA*sigma_c_ema + (1-ALPHA_EMA)*sigma_c_meas;
    sigma_c_ema      = max(sigma_c_ema, 0.8*sigma_c_platform);
    sigma_c_total_arr(f_idx) = sigma_c_ema;
    sigma_c_imu    = sigma_c_ema;
    sigma_c_final  = max(sigma_c_ema, 0.3);
    sigma_c_all(f_idx) = sigma_c_final;

    tp_n    = tgt(1).coor0 + t_fc*tgt(1).vt;
    los_vec = tp_n - ant_pos_fc;
    e_los_n = (los_vec/norm(los_vec)).';
    tgt_R_n = norm(los_vec);
    tgt_bin = min(max(round(tgt_R_n/range_res),1), N_fast);
    tgt_vr  = dot(tgt(1).vt, e_los_n.');
    [~, tgt_dop_guess] = min(abs(v_tf - (-tgt_vr)));

    % ---- A: MTI ----
    b_mti = [1, -1];
    frm_A_mti = filter(b_mti, 1, frm_raw_comp, [], 1);
    frm_A_mti(1, :) = 0;
    rd_A = fftshift(fft(frm_A_mti .* win_doppler, N_slow, 1), 1);

    % ---- B: Wide notch ----
    H_B  = make_wide_notch(v_tf, mu_c_final, WIDE_NOTCH_BW, N_slow, N_fast);
    rd_B = rd_pre .* H_B;

    % ---- C: RD adaptive soft notch ----
    H_C  = adaptive_soft_notch(v_tf, mu_c_final, B_flat_final, B_roll_final, ...
                               tfc.alpha, tfc.Hmin, N_slow, N_fast);
    rd_C0 = rd_pre .* H_C;

    cand_mask = detect_candidate_targets(abs(rd_C0).^2, 12);
    protect_cells = zeros(0,2);
    cl = find(cand_mask).';
    for k = cl
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

    % ---- D: TF ridge ----
    clut_lo = 50;
    clut_hi = max(clut_lo+1, tgt_bin-30);
    tfc.clut_bins = clut_lo:min(clut_hi, N_fast);

    B_core_f = max(BETA_CONSERVATIVE * B_flat_final, 0.5);
    rd_pass1 = apply_conservative_notch(rd_pre, v_tf, mu_c_final, B_core_f);

    [rd_D_pure, diagD] = tf_ridge_filter(frm_moco_comp, win_doppler, tfc, ...
                                          dv_att_pulse, sigma_c_imu, tgt_bin);

    cand_mask_d = detect_candidate_targets(abs(rd_pass1).^2,10);
    protect_cells_D = zeros(0,2);
    n_cand_raw_D  = 0;
    n_cand_kept_D = 0;

    cl = find(cand_mask_d).';
    B_flat_for_gate = max(B_flat_final, diagD.B_flat_used_ms);
    for k = cl
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
        dd = protect_cells_D(ci,1);
        rr = protect_cells_D(ci,2);
        dr = max(1,dd-2):min(N_slow,dd+2);
        rrg = max(1,rr-3):min(N_fast,rr+3);
        rd_D(dr, rrg) = rd_pre(dr, rrg);
    end

    % ---- Target peak ----
    [d_ref, r_ref, ~, ~, ok] = find_target_peak_prior(rd_pre, tgt_dop_guess, tgt_bin, ...
                                                       W_d_search, W_r_search, gamma_min_dB);
    peak_found(f_idx) = ok;
    if ~ok
        d_ref = tgt_dop_guess;
        r_ref = tgt_bin;
    end

    % ---- Metrics ----
    [SCR_A(f_idx), P_tgt_A]    = compute_scr(rd_A, d_ref, r_ref, N_slow, N_fast, clut_bins_est);
    [SCR_B(f_idx), P_tgt_B]    = compute_scr(rd_B, d_ref, r_ref, N_slow, N_fast, clut_bins_est);
    [SCR_C(f_idx), P_tgt_C]    = compute_scr(rd_C, d_ref, r_ref, N_slow, N_fast, clut_bins_est);
    [SCR_D(f_idx), P_tgt_D]    = compute_scr(rd_D, d_ref, r_ref, N_slow, N_fast, clut_bins_est);
    r_teri = max(1,r_ref-1):min(N_fast,r_ref+1);
    d_teri = max(1,d_ref-1):min(N_slow,d_ref+1);
    P_tgt_pre = mean(abs(rd_pre(d_teri,r_teri)).^2, 'all') + eps;
    TERI_A(f_idx) = P_tgt_A / P_tgt_pre;
    TERI_B(f_idx) = P_tgt_B / P_tgt_pre;
    TERI_C(f_idx) = P_tgt_C / P_tgt_pre;
    TERI_D(f_idx) = P_tgt_D / P_tgt_pre;
    tgt_separable(f_idx) = abs(-tgt_vr - mu_c_final) > B_flat_final + DELTA_MARGIN;

    BW3_C(f_idx) = compute_bw3dB(rd_C, tfc.clut_bins, v_tf, N_slow);
    BW3_D(f_idx) = compute_bw3dB(rd_D, tfc.clut_bins, v_tf, N_slow);

    % ---- CFAR & ROC ----
    [det_A, Z_A] = my_go_cfar(abs(rd_A).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_B, Z_B] = my_go_cfar(abs(rd_B).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_C, Z_C] = my_go_cfar(abs(rd_C).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    [det_D, Z_D] = my_go_cfar(abs(rd_D).^2, Tr_go,Gr_go,Td_go,Gd_go,Pfa_work);
    pts_A(f_idx) = sum(det_A,'all'); pts_B(f_idx) = sum(det_B,'all');
    pts_C(f_idx) = sum(det_C,'all'); pts_D(f_idx) = sum(det_D,'all');

    tgt_zone = false(N_slow, N_fast);
    tgt_zone(max(1,d_ref-2):min(N_slow,d_ref+2), ...
             max(1,r_ref-2):min(N_fast,r_ref+2)) = true;
    tgt_zone = tgt_zone & valid_cfar_mask;
    tgt_in_valid = any(tgt_zone, 'all');
    if tgt_in_valid
        if tgt_separable(f_idx)
            n_sep_frames = n_sep_frames + 1;
        else
            n_clut_frames = n_clut_frames + 1;
        end
    end

    Z_list = {Z_A, Z_B, Z_C, Z_D};
    pwr_list = {abs(rd_A).^2, abs(rd_B).^2, abs(rd_C).^2, abs(rd_D).^2};
    for s_idx = 1:4
        Z_s = Z_list{s_idx};
        p_s = pwr_list{s_idx};
        for p_idx = 1:nPfa
            alpha_p = N_ref_cfar*(PFA_SWEEP(p_idx)^(-1/N_ref_cfar)-1);
            det_p = (p_s > alpha_p*Z_s/N_ref_cfar) & valid_cfar_mask;
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
        end
    end
end

%% =========================================================================
% PART 10：计算最终指标
% =========================================================================
[~, pfa_ref_idx] = min(abs(PFA_SWEEP - Pfa_work));
Pd = pd_count(pfa_ref_idx, :) / num_frames;
FAR_A = fa_roc_cnt(pfa_ref_idx,1) / (num_frames * n_valid_cfar);
FAR_B = fa_roc_cnt(pfa_ref_idx,2) / (num_frames * n_valid_cfar);
FAR_C = fa_roc_cnt(pfa_ref_idx,3) / (num_frames * n_valid_cfar);
FAR_D = fa_roc_cnt(pfa_ref_idx,4) / (num_frames * n_valid_cfar);

Pd_roc  = pd_count / num_frames;
FAR_roc = fa_roc_cnt / (num_frames * n_valid_cfar);

n_pd_total_sep = n_sep_frames + n_clut_frames;
if n_sep_frames > 0
    Pd_sep_roc = pd_sep_cnt / n_sep_frames;
else
    Pd_sep_roc = zeros(nPfa, 4);
end
if n_clut_frames > 0
    Pd_clut_roc = pd_clut_cnt / n_clut_frames;
else
    Pd_clut_roc = zeros(nPfa, 4);
end

AUC = zeros(1, 4);
for s_idx = 1:4
    x = FAR_roc(:, s_idx);
    y = Pd_roc(:, s_idx);
    vld = isfinite(x) & isfinite(y) & x > 0;
    if sum(vld) >= 2
        AUC(s_idx) = trapz(log10(x(vld)), y(vld)) / ...
            (log10(max(x(vld))) - log10(min(x(vld))));
    end
end
vmask = peak_found;
nv = sum(vmask);
sep_mask = tgt_separable & vmask;

% ---- Pack results ----
results = struct();
results.SCR_mean = [mean(SCR_A), mean(SCR_B), mean(SCR_C), mean(SCR_D)];
if nv > 0
    results.SCR_mean_valid = [mean(SCR_A(vmask)), mean(SCR_B(vmask)), mean(SCR_C(vmask)), mean(SCR_D(vmask))];
    results.TERI = [mean(TERI_A(vmask)), mean(TERI_B(vmask)), mean(TERI_C(vmask)), mean(TERI_D(vmask))];
    if sum(sep_mask) > 5
        results.TERI_clean = [mean(TERI_A(sep_mask)), mean(TERI_B(sep_mask)), ...
                              mean(TERI_C(sep_mask)), mean(TERI_D(sep_mask))];
    else
        results.TERI_clean = nan(1,4);
    end
else
    results.SCR_mean_valid = nan(1,4);
    results.TERI = nan(1,4);
    results.TERI_clean = nan(1,4);
end
results.Pd = Pd;
results.FAR = [FAR_A, FAR_B, FAR_C, FAR_D];
results.AUC = AUC;
results.Pd_roc = Pd_roc;
results.FAR_roc = FAR_roc;
results.Pd_sep_roc = Pd_sep_roc;
results.Pd_clut_roc = Pd_clut_roc;
results.n_sep_frames = n_sep_frames;
results.n_clut_frames = n_clut_frames;
results.BW3 = [nan, nan, mean(BW3_C(BW3_C>0)), mean(BW3_D(BW3_D>0))];
results.SCR_gain = [0, mean(SCR_B)-mean(SCR_A), mean(SCR_C)-mean(SCR_A), mean(SCR_D)-mean(SCR_A)];

end  % end of main function

%% =========================================================================
%  ====================  内 置 函 数 区  ===================================
%% =========================================================================

% ─────────────────────────────────────────────────────────────────────────
function rd_out = apply_conservative_notch(rd_in, v_axis, mu_c, B_core)
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

k_ridge = viterbi_ridge(Ptf, k_pred, tfc.dk_search, tfc.sigmaT, tfc.lamT);

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

diag.Ptf=Ptf; diag.Ptf_tgt=Ptf_tgt; diag.k_pred=k_pred; diag.k_ridge=k_ridge; diag.G=G; diag.Nf=Nf;
diag.v_ridge = v_tf(k_ridge); diag.B_flat_v = B_flat_bin*dv_bin; diag.B_roll_v = B_roll_bin*dv_bin;
diag.w_r = w_r; diag.CNR_excess = CNR_excess;
diag.B_flat_used_ms = mean(B_flat_bin)*dv_bin;
diag.B_flat_max_ms  = max(B_flat_bin)*dv_bin;
end

% ─────────────────────────────────────────────────────────────────────────
function ridge = viterbi_ridge(Ptf, k_pred, dk, sigmaT, lamT)
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
H_vec = ones(N_dop,1);
H_vec(abs(v_axis(:)-center) <= half_bw) = 0;
H = repmat(H_vec, 1, N_rng);
end

% ─────────────────────────────────────────────────────────────────────────
function R = build_R_zyx_local(phi, theta, psi)
Rz = [cos(psi),-sin(psi),0; sin(psi),cos(psi),0; 0,0,1];
Ry = [cos(theta),0,sin(theta); 0,1,0; -sin(theta),0,cos(theta)];
Rx = [1,0,0; 0,cos(phi),-sin(phi); 0,sin(phi),cos(phi)];
R  = Rz*Ry*Rx;
end

% ─────────────────────────────────────────────────────────────────────────
function s_clut = sim_dynamic_clutter(Num_chirp, N_fast, ~, fc, PRI, K, Tc, ...
    t_fast, ant_pos_all, ant_vel_all, lambda, tgt_coor0, SCR_db)
% 注：调用前应在外部设置 rng 种子（MC 框架中固定为 42）
N_scatter = 80; c = 3e8; tgt_R0 = norm(tgt_coor0);
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
if nargin < 6 || isempty(clut_bins)
    clut_bins = 1:N_rng;
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
pwr = zeros(N_dop,1);
for b = clut_bins
    if b <= size(rd,2); pwr = pwr + abs(rd(:,b)).^2; end
end
pwr_n = pwr/(max(pwr)+eps);
bw = sum(pwr_n>=0.1)*(v_axis(2)-v_axis(1));
end

% ─────────────────────────────────────────────────────────────────────────
function [det, Z_map] = my_go_cfar(pwr_map, Tr, Gr, Td, Gd, Pfa)
[N_dop, N_rng] = size(pwr_map);
N_ref = 2*max(Td, Tr);   % GO-CFAR: ref cells in larger direction ×2
alpha_cfar = N_ref*(Pfa^(-1/N_ref)-1);

kernel_dop = [ones(Td,1); zeros(2*Gd+1,1); ones(Td,1)];
Z_d = conv2(pwr_map, kernel_dop, 'same');
kernel_rng = [ones(1,Tr), zeros(1,2*Gr+1), ones(1,Tr)];
Z_r = conv2(pwr_map, kernel_rng, 'same');
Z_map = max(Z_d, Z_r);
Z_map(Z_map < eps) = eps;

valid_mask = false(N_dop, N_rng);
dv_lo = Td+Gd+1; dv_hi = N_dop-Td-Gd;
rv_lo = Tr+Gr+1; rv_hi = N_rng-Tr-Gr;
valid_mask(dv_lo:dv_hi, rv_lo:rv_hi) = true;
Z_map(~valid_mask) = NaN;

det = (pwr_map > alpha_cfar*Z_map/N_ref) & valid_mask;
end

% ─────────────────────────────────────────────────────────────────────────
function [cand_mask, cand_list] = detect_candidate_targets(rd_power, thresh_db)
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
