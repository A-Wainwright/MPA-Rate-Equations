function plotDissociationResults(plot_data, params, sampleName)
    % Use sampleName in figure name
    figure('Name', [sampleName ' - Photodissociation Dynamics'], 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    
    t_fs = plot_data.t * 1e15;
    % To ensure graphs don't end up jagged, implement smoothing functions
    pulse_end_fs = plot_data.pulse_end_time * 1e15;
    window_size = max(5, round(length(t_fs) / 200));
    n0_smooth = movmean(plot_data.n0_t, window_size);
    n1_smooth = movmean(plot_data.n1_t, window_size);
    n2_smooth = movmean(plot_data.n2_t, window_size);   
    n0_recovered_smooth = movmean(plot_data.n0_1_t + plot_data.n0_2_t, window_size);
    dissociated_smooth = movmean(plot_data.dissociated_t, window_size);

    plot(t_fs, n0_smooth, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ground S_0 (smoothed)');
    hold on;
    plot(t_fs, n1_smooth, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Excited S_1');
    plot(t_fs, n2_smooth, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Excited S_2');
    plot(t_fs, n0_recovered_smooth, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ground (recovered)');
    plot(t_fs, dissociated_smooth, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Dissociated');
    % Mark end of pulse 
    xline(pulse_end_fs, 'r--', 'Alpha', 0.7, 'LineWidth', 2, 'DisplayName', 'Pulse End');
    xlabel('Time (fs)');
    ylabel('Population');
    title(sprintf('%s: Molecular Population Dynamics\nλ=%.0fnm, I=%.1eW/cm², τ=%.0ffs', ...
        sampleName, params.wavelength, params.I, params.tau_p));
    legend('show', 'Location', 'best');
    grid on;
    % Axis limits
    xlim([min(t_fs), max(t_fs)]);
    ylim([0, 1]);
end
