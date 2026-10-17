function manifest = prepare_real_imu(imu_dir, out_dir, max_files)
%PREPARE_REAL_IMU Convert recorded IMU CSV files to v4 replay CSV files.
%
% The source files contain, in order:
%   TID, ax, ay, az, gx, gy, gz, pitch, roll, yaw,
%   q0, q1, q2, q3, sample timestamp, sync timestamp
%
% The v4 experiment expects:
%   t, P_x, P_y, P_z, A_x(roll), A_y(pitch), A_z(yaw)
%
% Position is reconstructed as a nominal 10 m/s trajectory plus a
% detrended double integration of the measured acceleration. This keeps
% the recorded attitude dynamics while preventing long-term integration
% drift from dominating the replay.

if nargin < 1 || isempty(imu_dir)
    imu_dir = 'C:\Users\LiuShuo\Desktop\results\IMU';
end
if nargin < 2 || isempty(out_dir)
    this_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(this_dir, '..', 'data', 'imu_replay');
end
if nargin < 3 || isempty(max_files)
    max_files = inf;
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

files = dir(fullfile(imu_dir, '*.csv'));
if isempty(files)
    error('No IMU CSV files found in %s', imu_dir);
end

rows = cell(0, 1);
manifest = table();

for k = 1:min(numel(files), max_files)
    src = fullfile(files(k).folder, files(k).name);
    raw = readmatrix(src, 'NumHeaderLines', 1);
    raw = raw(~isnan(raw(:, 1)), :);
    if size(raw, 2) < 15
        warning('Skipping %s: fewer than 15 numeric columns.', files(k).name);
        continue;
    end

    tid = unwrap_tid(raw(:, 1));
    t = (tid - tid(1)) * 1e-3;
    if numel(t) < 2 || t(end) < 2.0
        warning('Skipping %s: replay duration is shorter than 2 s.', files(k).name);
        continue;
    end

    acc_body = raw(:, 2:4);
    pitch = raw(:, 8);
    roll  = raw(:, 9);
    yaw   = raw(:, 10);
    eul_deg = unwrap(deg2rad([roll, pitch, yaw]), [], 1) * 180 / pi;
    A_x = eul_deg(:, 1);
    A_y = eul_deg(:, 2);
    A_z = eul_deg(:, 3);

    acc_world = zeros(size(acc_body));
    g0 = 9.80665;
    for i = 1:size(acc_body, 1)
        R = rot_zyx(deg2rad(A_x(i)), deg2rad(A_y(i)), deg2rad(A_z(i)));
        acc_world(i, :) = (R * acc_body(i, :).') - [0; 0; g0];
    end

    % Remove slowly varying bias and smooth high-frequency integration noise.
    acc_world = detrend(acc_world, 'linear');
    win = max(3, round(0.05 / median(diff(t))));
    acc_world = movmean(acc_world, win, 1);
    v_res = cumtrapz(t, acc_world);
    v_res = detrend(v_res, 'linear');
    p_res = cumtrapz(t, v_res);
    p_res = detrend(p_res, 'linear');

    heading0 = deg2rad(A_z(1));
    v_nom = 10;
    P_x = v_nom * cos(heading0) * t + p_res(:, 1);
    P_y = v_nom * sin(heading0) * t + p_res(:, 2);
    P_z = p_res(:, 3);

    T = table(t, P_x, P_y, P_z, A_x, A_y, A_z);
    [~, stem, ~] = fileparts(files(k).name);
    out_file = fullfile(out_dir, ['Data_deduced_in_reverse_' stem '.csv']);
    writetable(T, out_file);

    rows{end+1, 1} = {files(k).name, out_file, t(end), ...
        std(A_x), std(A_y), std(A_z)}; %#ok<AGROW>
end

if ~isempty(rows)
    manifest = cell2table(vertcat(rows{:}), ...
        'VariableNames', {'source_file', 'replay_file', 'duration_s', ...
        'std_roll_deg', 'std_pitch_deg', 'std_yaw_deg'});
    manifest_file = fullfile(out_dir, 'imu_replay_manifest.csv');
    writetable(manifest, manifest_file);
    fprintf('[prepare_real_imu] Wrote %d replay files.\n', height(manifest));
    fprintf('[prepare_real_imu] Manifest: %s\n', manifest_file);
else
    warning('No replay files were generated.');
end
end

function tid = unwrap_tid(tid)
% The device counter wraps from 60000 to 1 in some files.
for i = 2:numel(tid)
    if tid(i) - tid(i-1) < -30000
        tid(i:end) = tid(i:end) + 60000;
    end
end
end

function R = rot_zyx(phi, theta, psi)
Rz = [cos(psi), -sin(psi), 0;
      sin(psi),  cos(psi), 0;
      0,         0,        1];
Ry = [cos(theta), 0, sin(theta);
      0,          1, 0;
      -sin(theta),0, cos(theta)];
Rx = [1, 0, 0;
      0, cos(phi), -sin(phi);
      0, sin(phi),  cos(phi)];
R = Rz * Ry * Rx;
end
