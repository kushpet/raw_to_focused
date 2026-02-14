# Range Doppler Algorithm (RDA) - Демонстрационный Набор

Реализация алгоритмов обработки SAR данных из **Cumming & Wong, Chapter 6**

## 📁 Файлы в наборе

### 1. **RDA_demo_cumming_ch6.m** - Основная демонстрация
Полная реализация трёх вариантов RDA из Figure 6.1:
- Basic RDA (рис. 6.1a)
- RDA with Accurate SRC - Option 2 (рис. 6.1b)  
- RDA with Approximate SRC - Option 3 (рис. 6.1c)

**Использует параметры из Table 6.1 книги Cumming (стр. 228)**

**Как запустить:**
```matlab
run RDA_demo_cumming_ch6.m
```
При запуске выберите вариант алгоритма (1, 2 или 3).

**Результат:**
- Ground truth положения целей
- Сырые SAR данные
- Визуализация каждого этапа обработки
- Сфокусированное изображение с метками истинных позиций

### 2. **RDA_single_target_demo.m** - Упрощённая демонстрация
Обработка ОДНОЙ точечной цели с детальной визуализацией каждого шага.

**Идеально для изучения физики алгоритма!**

**Как запустить:**
```matlab
run RDA_single_target_demo.m
```

**Что показывает:**
- **Этап 0**: Сырые данные (гиперболическая траектория цели)
- **Этап 1**: Range Compression (сжатие по дальности)
- **Этап 2**: Range-Doppler Domain (FFT по азимуту)
- **Этап 3**: RCMC (коррекция миграции)
- **Этап 4**: Azimuth Compression (фокусировка по азимуту)
- Анализ качества (PSLR по range и azimuth)

**Рекомендуется начать с этого файла!**

### 3. **RDA_physics_tutorial.m** - Физика и теория
Подробные объяснения всех концепций на русском языке:

**Содержание:**
- Что такое Doppler Centroid и почему он важен
- Влияние Squint Angle на f_dc
- Зачем делают FFT по азимуту
- Физический смысл Range-Doppler координат
- Что такое Range Cell Migration и как его корректировать
- Зачем нужен Secondary Range Compression (SRC)
- Пошаговый разбор всех вариантов алгоритма
- Практические замечания
- Связь с алгоритмами Sentinel-1

**Как использовать:**
Откройте файл в редакторе MATLAB и читайте секции.  
Это НЕ исполняемый код, а справочная документация.

### 4. **RDA_formulas_reference.m** - Математический справочник
Сборник всех ключевых формул из Cumming Chapter 6 с указанием страниц и номеров уравнений.

**Содержание:**
- Геометрические соотношения
- Формулы Range Compression
- Range-Doppler domain
- Range Cell Migration
- RCMC операции
- Azimuth Compression
- Secondary Range Compression (оба варианта)
- Импульсная характеристика
- Численные примеры из Table 6.1
- Критерии применимости алгоритмов

**Как использовать:**
Справочник для быстрого поиска формул и констант.

## 🎯 Рекомендуемый порядок изучения

1. **Сначала прочитайте теорию:**
   ```matlab
   edit RDA_physics_tutorial.m
   ```
   Разделы 1.1-1.5 дают понимание ключевых концепций.

2. **Запустите упрощённую демонстрацию:**
   ```matlab
   run RDA_single_target_demo.m
   ```
   Смотрите на картинки, понимайте, что происходит на каждом этапе.

3. **Поэкспериментируйте с параметрами:**
   Откройте `RDA_single_target_demo.m`, измените:
   - `theta_squint` → увидите изменение f_dc и траектории
   - `R0` → изменит величину миграции
   - `Vr` → повлияет на Doppler и миграцию
   - `SNR_dB = 10` → увидите влияние шума

4. **Изучите полную реализацию:**
   ```matlab
   run RDA_demo_cumming_ch6.m
   ```
   Выберите вариант алгоритма. Сравните результаты разных вариантов.

5. **Справочники всегда под рукой:**
   - `RDA_physics_tutorial.m` - для понимания физики
   - `RDA_formulas_reference.m` - для формул

## 🔑 Ключевые концепции (краткая справка)

### Doppler Centroid (f_dc)
```
f_dc = -(2·V/λ)·sin(θ_squint)
```
- При squint = 0° (broadside): f_dc = 0 Гц
- При squint = 3.5°: f_dc ≈ 320 Гц
- Если не компенсировать → расфокусировка по азимуту

### Range Cell Migration (RCM)
```
ΔR(f_η) = R0·[sqrt(1 - (λ·f_η/(2·V))²) - 1]
```
- Цель движется по гиперболе → сигнал мигрирует через range bins
- В Range-Doppler domain: каждая f_η имеет свою миграцию
- RCMC исправляет это интерполяцией

### Secondary Range Compression (SRC)
```
H_SRC(f_r, f_η) = exp(j·π·f_r²/Km(f_η))
где Km(f_η) = -λ²·R0/(8·V²·D(f_η)³)
```
- Корректирует остаточную coupling между range и azimuth
- **Option 2**: Точный, применяется ПЕРЕД azimuth compression
- **Option 3**: Приближённый, применяется ПОСЛЕ azimuth compression

## 📊 Что смотреть на визуализациях

### Raw Data (сырые данные)
- **Magnitude**: Гиперболические траектории целей
- **Phase**: Видна chirp модуляция
- Энергия размазана по range и azimuth

### Range Compressed
- **Magnitude**: Траектории стали уже по range
- Range resolution улучшилось до ρ_r = c/(2·Br)
- Азимутальная расфокусировка осталась

### Range-Doppler Domain
- **Magnitude**: Энергия целей на разных Doppler частотах
- Видна миграция по range для разных f_η
- Это ключевой домен для RCMC!

### After RCMC
- **Magnitude**: Траектории выровнены (горизонтальные линии)
- Энергия цели собрана в один range bin
- Готово для azimuth compression

### Focused Image
- **Magnitude**: Точечные цели сфокусированы
- Оцениваем качество: PSLR, разрешение
- Сравниваем с ground truth

## 🔧 Требования

- **MATLAB R2019b** или новее
- **Toolboxes**: 
  - Signal Processing Toolbox (для FFT/IFFT)
  - Phased Array System Toolbox (только для полной демо с симуляцией)

Для упрощённой демонстрации (`RDA_single_target_demo.m`) достаточно базового MATLAB.

## 📖 Ссылки на литературу

### Основные источники

**[1] Cumming I.G., Wong F.H.**  
"Digital Processing of Synthetic Aperture Radar Data: Algorithms and Implementation"  
Artech House, 2005  
**Chapter 6: The Range Doppler Algorithm** (pages 225-262)

Ключевые разделы:
- Section 6.1: Introduction (стр. 225)
- Section 6.2: Algorithm Overview (стр. 226)
- Section 6.3: RDA in the Low Squint Case (стр. 229)
  - 6.3.1: Range Compression (стр. 230)
  - 6.3.2: Azimuth FFT (стр. 231)
  - 6.3.3: RCMC (стр. 233)
  - 6.3.4: Azimuth Compression (стр. 236)
- Section 6.4: Secondary Range Compression (стр. 244)
  - 6.4.2: Option 2 - Before Azimuth Compression (стр. 248)
  - 6.4.3: Option 3 - After Azimuth Compression (стр. 250)

**[2] Sentinel-1 Level 1 Detailed Algorithm Definition**  
Document: S1-TN-MDA-52-7445, Issue 2.5, 2022  
https://sentinel.esa.int/documents/247904/1653442/Sentinel-1-Level-1-Detailed-Algorithm-Definition

Релевантные разделы:
- Section 4: SAR Focusing
- Section 4.2: Range Doppler Algorithm
- Section 4.3: Extended Chirp Scaling

### Дополнительная литература

**[3] Cumming I.G., Neo Y.L., Wong F.H.**  
"Interpretations of the Omega-K algorithm and comparisons with other algorithms"  
IGARSS 2003

**[4] Davidson G.W., Cumming I.G.**  
"Signal Properties of Spaceborne Squint-mode SAR"  
IEEE TGRS, Vol. 35, No. 3, May 1997

**[5] Raney R.K., Runge H., Bamler R., Cumming I.G., Wong F.H.**  
"Precision SAR processing using chirp scaling"  
IEEE TGRS, Vol. 32, No. 4, July 1994

## 🎓 Следующие шаги

После освоения базового RDA:

1. **Изучите модификации для больших squint углов**
   - Extended Chirp Scaling (ECS)
   - Omega-K algorithm

2. **Разберитесь с Sentinel-1 IPF алгоритмами**
   - Preprocessing: Doppler centroid estimation
   - Timing и orbit corrections
   - Radiometric calibration

3. **Перейдите к реальным данным**
   - Загрузите Sentinel-1 Level-0 данные
   - Примените изученные алгоритмы
   - Сравните с официальными Level-1 продуктами

## 🐛 Известные ограничения

- Упрощённая симуляция (не учитывает):
  - Rotation Earth effects
  - Orbit variations
  - Atmospheric effects
  - Antenna pattern реалистичный
  
- Basic RDA работает лучше при:
  - Малый squint (< 5°)
  - Узкий swath
  - Средняя разрешение

- Для реальных Sentinel-1 данных нужны дополнительные этапы:
  - Orbit state vectors interpolation
  - Precise timing corrections
  - Interferometric registration

## 💡 Советы по отладке

Если изображение расфокусировано:

1. **Проверьте Doppler Centroid:**
   ```matlab
   fprintf('f_dc = %.1f Hz\n', f_dc);
   ```
   Должен соответствовать геометрии.

2. **Проверьте RCMC:**
   После RCMC траектория должна быть горизонтальной.
   Посмотрите на `stages.rcmc` визуализацию.

3. **Проверьте SRC:**
   Если цели на разных дистанциях фокусируются по-разному,
   нужен SRC (алгоритм 2 или 3).

4. **Проверьте параметры:**
   ```matlab
   fprintf('Max migration: %.3f м\n', max_migration);
   fprintf('Range resolution: %.3f м\n', rho_r);
   ```
   Если миграция > ρ_r/2 → RCMC критичен.

## 📧 Обратная связь

Это демонстрационный набор для изучения RDA.  
Не предназначен для обработки production данных.

Для работы с реальными Sentinel-1 данными см.:
- SNAP (Sentinel Application Platform)
- ISCE (InSAR Scientific Computing Environment)
- SARscape

---

**Удачи в изучении SAR обработки! 🛰️📡**
