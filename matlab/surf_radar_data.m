function surf_radar_data(Z, decim)
% SURF_RADAR_DATA  Отображает амплитуду комплексной матрицы Z с прореживанием.
%   surf_radar_data(Z, decim)
%   Z     — комплексная (или вещественная) матрица
%   decim — шаг прореживания по обоим измерениям (целое >= 1)

    if nargin < 2 || isempty(decim)
        decim = 1;
    end
    if decim < 1 || decim ~= floor(decim)
        error('decim должен быть целым числом >= 1');
    end

    % Амплитуда
    A = abs(Z);

    % Прореживание
    Adec = A(1:decim:end, 1:decim:end);

    % Координаты (можно заменить на физические, если известны)
    [nx, ny] = size(Adec);
    [X, Y] = meshgrid(1:ny, 1:nx);

    % Отображение
    figure;
    h = surf(X, Y, Adec);      % [web:2]
    set(h, 'EdgeColor', 'none');  % убрать рёбра [web:2][web:4]
    colormap('jet');           % цветная картинка
    colorbar;

    xlabel('X');
    ylabel('Y');
    zlabel('Amplitude');
    title('Radar amplitude surf');

    view(3);
end

