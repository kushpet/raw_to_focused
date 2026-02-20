% =========================================================================
% Simplified RDA Demo - Single Point Target
% =========================================================================
%
% Упрощённая демонстрация для изучения каждого этапа RDA
% на примере ОДНОЙ точечной цели
%
% Позволяет:
%   - Визуализировать сигнал на каждом этапе
%   - Понять влияние каждой операции
%   - Экспериментировать с параметрами
%
% =========================================================================

clear; close all; clc;

%% ========================================================================
% ПАРАМЕТРЫ - упрощённые для наглядности
% =========================================================================

c = 3e8;                            % Скорость света [м/с]
f0 = 5.3e9;                         % Несущая частота [Гц]
lambda = c / f0;                    % Длина волны [м]

% Параметры платформы
Vr = 150;                           % Скорость [м/с]

% Параметры сигнала
Tr = 2.5e-6;                        % Длительность импульса [с]
Kr = 20e12;                         % Частота модуляции [Гц/с]
Br = abs(Kr * Tr);                  % Полоса [Гц]

% Параметры дискретизации
Fs = 60e6;                          % Частота дискретизации по range [Гц]
PRF = 1200;                         % ЧПИ [Гц]

% Параметры сцены
R0 = 2e3;                           % Дистанция до цели [м]
theta_squint = deg2rad(3.5);        % Угол squint [рад]

% Размеры данных (малые для наглядности)
Naz = 2560;                         % Число импульсов
Nrg = 320;                          % Отсчётов на импульс

Lsca = Naz / PRF * Vr;              % Размер сцены по азимуту [м]
La = lambda * R0 / Lsca;            % Апертура антенны в азимутуальной плоскости [м]

% Доплеровские параметры
f_dc = -2*Vr/lambda * sin(theta_squint); % Допплер центроид [Гц]

fprintf('=== ПАРАМЕТРЫ ===\n');
fprintf('Несущая частота: %.2f ГГц\n', f0/1e9);
fprintf('Длина волны: %.4f м\n', lambda);
fprintf('Разрешение по дальности: %.2f м\n', c/(2*Br));
fprintf('Угол squint: %.2f град\n', rad2deg(theta_squint));
fprintf('Допплер центроид: %.1f Гц\n', f_dc);
fprintf('\n');

%% ========================================================================
% ГЕНЕРАЦИЯ СЫРЫХ ДАННЫХ - одна цель в центре
% =========================================================================

fprintf('Генерация сырых данных для одной цели...\n');

% Оси времени
tau0   = 2*R0/c;                               % задержка до цели
t_fast = tau0 + ((0:Nrg-1) - Nrg/2) / Fs;      % окно вокруг tau0 (строб приёма)
eta_slow = ((0:Naz-1) - Naz/2) / PRF; % Медленное время [с]

% Положение цели - в центре
eta0 = 0;                           % Цель в центре окна по медленному времени
A0 = 1.0;                           % Амплитуда

% Инициализация
S_raw = zeros(Nrg, Naz);

% Генерация сигнала
for m = 1:Naz
    eta = eta_slow(m);  % Время
    
    % История дальности (гипербола)
    R_eta = sqrt(R0^2 + Vr^2 * (eta - eta0)^2);
    
    % Задержка сигнала до цели и обратно
    tau = 2 * R_eta / c;
    
    % Угол на цель
    theta_target = atan2(Vr * (eta-eta0), R0);
    % Доплеровское смещение (справочно, не используется в расчётах)
    doppler = 2 * Vr / lambda * sin(theta_target);
    % Угол относительно максимума антенны
    theta_wa = theta_target - theta_squint;
    % ЛЧМ сигнал
    % по формуле 4.39 из Cumming, 2005.
    % t_fast - окно, фиксировано для R0
    % tau - задержка до цели и обратно, зависит от текущей дальности
    % t_centered - ноль соответствует середине LFM сигнала
    t_centered = t_fast - tau;
    % окно для вырезки LFM - аналог w_r
    range_window = double(abs(t_centered) <= Tr/2);
    % Учёт диаграммы направленности антенны
    w_a = sinc(La * sin(theta_wa) / lambda);   
    % Фазовый член - изменение фазы несущей
    phase_term = exp(-1j * 4 * pi * R_eta / lambda);
    % ЛЧМ по полному окну
    chirp = exp(1j * pi * Kr * t_centered.^2);
    
    % Сигнал
    S_raw(:, m) = A0 * w_a * phase_term * (range_window .* chirp).';
end

fprintf('Max |S_raw| = %.3f\n', max(abs(S_raw(:))));

% Добавить шум
SNR_dB = 25;
signal_power = sum(abs(S_raw(:)).^2) / numel(S_raw);
noise_power = signal_power / (10^(SNR_dB/10));
S_raw = S_raw + sqrt(noise_power/2) * (randn(size(S_raw)) + 1j*randn(size(S_raw)));

% Обнулить края для последующей циклической свёртки
Lref = Tr * Fs + 1;    % Длинна опорного сигнала (нечетная)
Blank_raw = zeros(floor(Lref/2), Naz);
S_raw(1:floor(Lref/2), :) = Blank_raw; % Начало по дальности
S_raw(Nrg - floor(Lref/2) + 1:Nrg, :) = Blank_raw; % Конец по дальности
% Здесь ожидаем максимум свёрнутого сигнала
center_range_bin = floor(Nrg/2) + 1;
center_azimuth_bin = floor(Naz/2) + 1;

fprintf('Сырые данные: %d × %d\n\n', Nrg, Naz);

%% ========================================================================
% ВИЗУАЛИЗАЦИЯ: Raw Data
% =========================================================================

figure('Name', 'ЭТАП 0: Сырые данные', 'Position', [50 50 1400 500]);
% Реальная часть
subplot(2,3,1);
imagesc(eta_slow*1e3, t_fast*1e6, real(S_raw));
xlabel('Медленное время [мс]'); ylabel('Быстрое время [мкс]');
title('Реальная часть сырых данных');
colorbar; axis xy; colormap(gca, jet);

% Мнимая часть
subplot(2,3,2);
imagesc(eta_slow*1e3, t_fast*1e6, imag(S_raw));
xlabel('Медленное время [мс]'); ylabel('Быстрое время [мкс]');
title('Мнимая часть сырых данных');
colorbar; axis xy; colormap(gca, jet);

% Амплитуда
subplot(2,3,4);
imagesc(eta_slow*1e3, t_fast*1e6, abs(S_raw));
xlabel('Медленное время [мс]'); ylabel('Быстрое время [мкс]');
title('Амплитуда сырых данных');
colorbar; axis xy; colormap(gca, jet);

% Фаза
subplot(2,3,5);
imagesc(eta_slow*1e3, t_fast*1e6, angle(S_raw));
xlabel('Медленное время [мс]'); ylabel('Быстрое время [мкс]');
title('Фаза сырых данных');
colorbar; axis xy; colormap(gca, hsv);

% Азимутальный профиль (центральный range bin)
subplot(2,3,6);
plot(eta_slow*1e3, abs(S_raw(center_range_bin, :)), 'b-', 'LineWidth', 1.5);
xlabel('Медленное время [мс]'); ylabel('Амплитуда');
title(sprintf('Азимутальный профиль (range bin %d)', center_range_bin));
grid on;

% Пояснение
annotation('textbox', [0.02 0.95 0.3 0.04], 'String', ...
    'Видна диаграмма направленности по азимуту', ...
    'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');

pause(1);

%% ========================================================================
% ЭТАП 1: RANGE COMPRESSION
% =========================================================================

fprintf('ЭТАП 1: Range Compression...\n');

% Опорный чирп от отрицательного до положетельного времени
etalon_chirp_time = ((0:Tr*Fs) - Tr*Fs/2) / Fs;
etalon_chirp = exp(1j * pi * Kr * etalon_chirp_time.^2);

% Нулепаддинг до Nrg
etalon_padded = zeros(Nrg, 1);
etalon_padded(1:min(Lref, Nrg)) = etalon_chirp;
% Сдвигаем циклически назад, так, чтобы центр LFM попал 
% в нулевую (первую) позицию. Тогда при выполнении
% взаимной корреляции получим пик в точке середины LFM 
% принятого сигнала
etalon_padded = circshift(etalon_padded, - floor(Lref / 2));
% Спектр опорного
etalon_spectr = fft(etalon_padded, Nrg);

% Циклическая свёртка (размер Nrg)
S_rc = zeros(Nrg, Naz);
for m = 1:Naz
    S_fft = fft(S_raw(:, m), Nrg);
    S_compressed = S_fft .* conj(etalon_spectr);
    S_rc(:, m) = ifft(S_compressed);
end

fprintf(' Разрешение по дальности: %.2f м\n\n', c/(2*Br));

%% ВИЗУАЛИЗАЦИЯ: Range Compressed
figure('Name', 'ЭТАП 1: Range Compression', 'Position', [60 60 1400 500]);

subplot(1,3,1);
imagesc(eta_slow*1e3, t_fast*1e6, abs(S_rc));
xlabel('Медленное время [мс]'); ylabel('Быстрое время [мкс]');
title('После Range Compression (пик в центр)');
colorbar; axis xy; colormap(gca, jet);

subplot(1,3,2);
plot(t_fast*1e6, abs(S_rc(:, center_azimuth_bin)), 'r-', 'LineWidth', 1.5);
xlabel('Быстрое время [мкс]'); ylabel('Амплитуда');
title('Range профиль (центральный импульс)');
grid on;

subplot(1,3,3);
plot(eta_slow*1e3, abs(S_rc(center_range_bin, :)), 'b-', 'LineWidth', 1.5);
xlabel('Медленное время [мс]'); ylabel('Амплитуда');
title(sprintf('Азимутальный профиль (range bin %d)', center_range_bin));
grid on;

annotation('textbox', [0.02 0.95 0.4 0.04], 'String', ...
    'Сжатие по дальности: цель локализована по range', ...
    'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');

pause(1);

%% ========================================================================
% ЭТАП 2: AZIMUTH FFT → Range-Doppler Domain
% =========================================================================

fprintf('ЭТАП 2: Переход в Range-Doppler домен...\n');

% FFT по азимуту
S_rd = fftshift(fft(S_rc, [], 2), 2);

% Ось Допплеровских частот
f_eta = ((0:Naz-1) - Naz/2) * (PRF/Naz) + f_dc;

fprintf('  Диапазон Допплеровских частот: [%.1f, %.1f] Гц\n\n', ...
        min(f_eta), max(f_eta));

%% ВИЗУАЛИЗАЦИЯ: Range-Doppler Domain
figure('Name', 'ЭТАП 2: Range-Doppler Domain', 'Position', [70 70 1400 500]);

subplot(1,3,1);
imagesc(f_eta, t_fast*1e6, abs(S_rd));
xlabel('Допплеровская частота [Гц]'); ylabel('Быстрое время [мкс]');
title('Range-Doppler Domain');
colorbar; axis xy; colormap(gca, jet);

subplot(1,3,2);
plot(f_eta, abs(S_rd(center_range_bin, :)), 'g-', 'LineWidth', 1.5);
xlabel('Допплеровская частота [Гц]'); ylabel('Амплитуда');
title('Допплеровский спектр (центральный range bin)');
grid on;

% Показать траекторию миграции
subplot(1,3,3);
imagesc(f_eta, t_fast*1e6, abs(S_rd));
xlabel('Допплеровская частота [Гц]'); ylabel('Быстрое время [мкс]');
title('Траектория Range Cell Migration');
colorbar; axis xy; colormap(gca, jet);
hold on;

% Теоретическая траектория миграции
for k = 1:length(f_eta)
    f_d = f_eta(k); % Просто частота, которая соответствует этому "бину"
    % Истинная дальность гиперболы для этого f_d
    sin_theta = lambda * f_d / (2 * Vr);  % sin(θ) = λ f_d / (2 V_r)
    cos_theta = sqrt(1 - sin_theta^2);    % cos(θ)
    R_actual = R0 / cos_theta;            % R(θ) = R0 / cos(θ) > R0
    delta_R = R_actual - R0;              % ПОЛОЖИТЕЛЬНОЕ!
    t_migration = 2 * R_actual / c;       % > tau0
    plot(f_eta(k), t_migration*1e6, 'r.', 'MarkerSize', 8);
end
hold off;
legend('Данные', 'Теоретическая траектория', 'Location', 'best');

annotation('textbox', [0.02 0.95 0.5 0.04], 'String', ...
    'Видна миграция цели по range bins для разных Допплер частот', ...
    'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');

pause(1);

%% ========================================================================
% ЭТАП 3: RANGE CELL MIGRATION CORRECTION (RCMC)
% =========================================================================

fprintf('ЭТАП 3: Range Cell Migration Correction...\n');

S_rcmc = zeros(Nrg, Naz);
for m = 1:Naz
    f_eta_m = f_eta(m);
    
    % 1. УГОЛ (точно)
    sin_theta = lambda * f_eta_m / (2 * Vr);
    cos_theta = sqrt(max(0, 1 - sin_theta^2));  % защита от |sin| > 1
    
    % 2. ДАЛЬНОСТЬ гиперболы (Cumming eq. 4.2-12)
    R_slant = R0 / cos_theta;  % > R0!
    
    % 3. Миграция (ПОЛОЖИТЕЛЬНАЯ)
    delta_R = R_slant - R0;
    delta_t = 2 * delta_R / c;
    
    % 4. ИНТЕРПОЛЯЦИЯ НАЗАД к R0
    t_shifted = t_fast + delta_t;  % Плюс!
    
    S_rcmc(:, m) = interp1(t_fast, S_rd(:, m), t_shifted, 'linear', 0);
end

fprintf('  Максимальная миграция: %.3f м\n\n', abs(min(R0 * (sqrt(1 - (lambda * f_eta / (2*Vr)).^2) - 1))));

%% ВИЗУАЛИЗАЦИЯ: After RCMC
figure('Name', 'ЭТАП 3: RCMC', 'Position', [80 80 1400 500]);

subplot(1,3,1);
imagesc(f_eta, t_fast*1e6, abs(S_rcmc));
xlabel('Допплеровская частота [Гц]'); ylabel('Быстрое время [мкс]');
title('После RCMC');
colorbar; axis xy; colormap(gca, jet);

subplot(1,3,2);
imagesc(f_eta, t_fast*1e6, abs(S_rd));
xlabel('Допплеровская частота [Гц]'); ylabel('Быстрое время [мкс]');
title('До RCMC (для сравнения)');
colorbar; axis xy; colormap(gca, jet);

subplot(1,3,3);
plot(f_eta, abs(S_rcmc(center_range_bin, :)), 'b-', 'LineWidth', 1.5);
xlabel('Допплеровская частота [Гц]'); ylabel('Амплитуда');
title('Допплеровский спектр после RCMC');
grid on;

annotation('textbox', [0.02 0.95 0.6 0.04], 'String', ...
    'Миграция исправлена: энергия цели выровнена по range (горизонтальная линия)', ...
    'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');

pause(1);

%% ========================================================================
% ЭТАП 4: AZIMUTH COMPRESSION
% =========================================================================

fprintf('ЭТАП 4: Azimuth Compression...\n');

% Азимутальный согласованный фильтр
H_ac = zeros(1, Naz);
for m = 1:Naz
    f_eta_m = f_eta(m);
    D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
    H_ac(m) = exp(1j * 4 * pi * R0 * D_f / lambda);
end

% Применить фильтр
S_ac = S_rcmc .* repmat(H_ac, Nrg, 1);

% Обратное FFT по азимуту
S_focused = ifft(ifftshift(S_ac, 2), [], 2);

fprintf('  Изображение сфокусировано!\n\n');

%% ВИЗУАЛИЗАЦИЯ: Focused Image
figure('Name', 'ЭТАП 4: Сфокусированное изображение', 'Position', [90 90 1400 500]);

subplot(1,3,1);
imagesc(eta_slow*1e3, t_fast*1e6, abs(S_focused));
xlabel('Азимут [мс]'); ylabel('Дальность [мкс]');
title('Сфокусированное изображение');
colorbar; axis xy; colormap(gca, jet);
hold on;
% Отметить истинное положение цели
plot(0, (2*R0/c)*1e6, 's-', 'MarkerSize', 20, 'LineWidth', 2, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
hold off;

% Range профиль
subplot(1,3,2);
plot(t_fast*1e6, abs(S_focused(:, Naz/2)), 'r-', 'LineWidth', 1.5);
xlabel('Дальность [мкс]'); ylabel('Амплитуда');
title('Range профиль через цель');
grid on;
xlim([t_fast(1)*1e6, t_fast(end)*1e6]);

% Азимутальный профиль
subplot(1,3,3);
plot(eta_slow*1e3, abs(S_focused(center_range_bin, :)), 'b-', 'LineWidth', 1.5);
xlabel('Азимут [мс]'); ylabel('Амплитуда');
title('Азимутальный профиль через цель');
grid on;
xlim([eta_slow(1)*1e3, eta_slow(end)*1e3]);

annotation('textbox', [0.02 0.95 0.4 0.04], 'String', ...
    'Цель полностью сфокусирована по range и азимуту', ...
    'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');

%% ========================================================================
% СРАВНЕНИЕ: До и После
% =========================================================================

figure('Name', 'СРАВНЕНИЕ: До и После фокусировки', 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(eta_slow*1e3, t_fast*1e6, abs(S_raw));
xlabel('Азимут [мс]'); ylabel('Дальность [мкс]');
title('ДО: Сырые данные');
colorbar; axis xy; colormap(gca, jet);

subplot(1,2,2);
imagesc(eta_slow*1e3, t_fast*1e6, abs(S_focused));
xlabel('Азимут [мс]'); ylabel('Дальность [мкс]');
title('ПОСЛЕ: Сфокусированное изображение');
colorbar; axis xy; colormap(gca, jet);
hold on;
plot(0, (2*R0/c)*1e6, 's-', 'MarkerSize', 20, 'LineWidth', 2, ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
hold off;

%% ========================================================================
% АНАЛИЗ КАЧЕСТВА ФОКУСИРОВКИ
% =========================================================================

fprintf('=== АНАЛИЗ КАЧЕСТВА ===\n');

% Найти пик
[max_val, max_idx] = max(abs(S_focused(:)));
[peak_range, peak_az] = ind2sub(size(S_focused), max_idx);

fprintf('Пик сигнала:\n');
fprintf('  Range bin: %d (ожидалось ~%d)\n', peak_range, center_range_bin);
fprintf('  Azimuth sample: %d (ожидалось ~%d)\n', peak_az, center_azimuth_bin);
fprintf('  Амплитуда: %.2f\n', max_val);

% PSLR (Peak to Side Lobe Ratio) по range
range_profile = abs(S_focused(:, peak_az));
range_profile = range_profile / max(range_profile); % Нормализация
[~, peak_idx] = max(range_profile);

% Исключить главный лепесток (±5 отсчётов)
mainlobe_width = 5;
sidelobe_region = range_profile;
sidelobe_region(max(1, peak_idx-mainlobe_width):min(Nrg, peak_idx+mainlobe_width)) = 0;
max_sidelobe = max(sidelobe_region);

PSLR_range = 20*log10(1 / max_sidelobe);
fprintf('\nКачество по Range:\n');
fprintf('  PSLR: %.2f дБ\n', PSLR_range);

% PSLR по азимуту
az_profile = abs(S_focused(peak_range, :));
az_profile = az_profile / max(az_profile);
[~, peak_idx_az] = max(az_profile);

sidelobe_region_az = az_profile;
sidelobe_region_az(max(1, peak_idx_az-mainlobe_width):min(Naz, peak_idx_az+mainlobe_width)) = 0;
max_sidelobe_az = max(sidelobe_region_az);

PSLR_az = 20*log10(1 / max_sidelobe_az);
fprintf('\nКачество по Азимуту:\n');
fprintf('  PSLR: %.2f дБ\n', PSLR_az);

%% ========================================================================
% ЭКСПЕРИМЕНТЫ
% =========================================================================

fprintf('\n=== ЭКСПЕРИМЕНТЫ ===\n');
fprintf('Попробуйте изменить параметры в начале файла:\n');
fprintf('  1. theta_squint - увидите изменение f_dc и траектории\n');
fprintf('  2. R0 - изменит величину миграции\n');
fprintf('  3. Vr - повлияет на Допплер и миграцию\n');
fprintf('  4. SNR_dB - уменьшите до 10 дБ и посмотрите на шум\n');
fprintf('\n');

fprintf('=== ГОТОВО ===\n');
fprintf('Все этапы RDA визуализированы.\n');
fprintf('Изучите картинки для понимания каждого шага.\n');
