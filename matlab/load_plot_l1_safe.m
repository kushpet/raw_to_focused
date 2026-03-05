%% Путь к SAFE-папке
%safeDir = 'C:\Users\Dmitry\Documents\GitHub\raw_to_focused\python\data\la_ocean\SLC\S1A_IW_SLC__1SDV_20260226T140055_20260226T140122_063391_07F66A_5FE2.SAFE';
%safeDir = 'C:\Users\Dmitry\Documents\GitHub\raw_to_focused\python\data\rosemond\SLC\S1A_IW_SLC__1SDV_20260209T135208_20260209T135235_063143_07ED08_7039.SAFE';
safeDir = 'D:\SarData\2026.02.18_Australia_Etalon\S1C_IW_SLC__1SSV_20260218T083150_20260218T083221_006407_00CE5B_2A1B.SAFE';
%% Разбор manifest.safe
maniHead = parseSafeManifest(safeDir);   % из sent1-L1-utilities [web:39]

% Индексы: IW1, VV
sidx = 1;   % subswath index: 1=IW1, 2=IW2, 3=IW3 (для IW продукта) [web:6]
pidx = 1;   % polarization index: 1=VV (первый в списке для DV‑продукта) [web:3]
vidx = 1;   % первый slice (обычно 1)

% Файл с данными SLC для IW1/VV
dataFile = fullfile(maniHead.rootDir, maniHead.dataFile{sidx,pidx,vidx});
fprintf('Using data file:\n%s\n', dataFile);

%% Чтение комплексных SLC-данных
% true = читать как комплексный массив (см. документацию readSent1Data) [web:39]
Z = readSent1Data(dataFile, true);   % complex single, размер [azimuth x range]

%% Амплитуда / фаза
A    = abs(Z);
A_dB = 20*log10(A + eps);
Ph   = angle(Z);

dsAz = 5;    % шаг по азимуту (строкам)
dsRg = 5;    % шаг по дальности (столбцам)

A_dec   = A(1:dsAz:end, 1:dsRg:end);
A_dec_dB = 20*log10(A_dec + 1.0);

%% 2D отображение амплитуды (dB)
figure('Color','w');
imagesc(A_dec);
axis image; colorbar;
xlabel('Range (samples)');
ylabel('Azimuth (lines)');
title('Sentinel-1 SLC IW1 VV amplitude (dB) (downsampled)');

%% При желании — фаза
% figure('Color','w');
% imagesc(Ph); axis image; colormap(hsv); colorbar;
% title('Sentinel-1 SLC IW1 VV phase (rad)');

%% 3D поверхность (амплитуда как высота, с даунсэмплингом)

figure('Color','w');
surf(A_dec, 'EdgeColor','none');
colormap(turbo); colorbar;
view(45,60);
xlabel('Range (samples)'); ylabel('Azimuth (lines)'); zlabel('Amplitude (dB)');
title('Sentinel-1 SLC IW1 VV amplitude surface (downsampled)');
%camlight headlight; lighting gouraud;
