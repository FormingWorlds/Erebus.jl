function [result_1, result_2] = compare_mat(func, data, range, text)
%compare_mat compare two data series wrt to func over all steps in a range
    result_1 = zeros(length(range));
    result_2 = zeros(length(range));
    for step = range
        [data_1, data_2] = load_mat(data, step);
        result_1(step) = func(data_1);
        result_2(step) = func(data_2);
    end
    h = figure;
    plot(range, result_1, '-*')
    hold on
    plot(range, result_2, '-o')
    hold off
    legend('MATLAB reference', 'Julia HydrologyPlanetesimals')
    title([data ': ' text])
    xlabel('time step')
    ylabel([get_var_name(func) '(' data ')'])
    saveas(h, [data '_' get_var_name(func) '.pdf']);
%     close(h);
end