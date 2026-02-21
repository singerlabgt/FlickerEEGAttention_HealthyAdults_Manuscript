function [] = plotfooof(fooofresult)
    %plotfooof Written by Matty 2024-12-09
    %   Adds subject


    % Plot FOOOF results:
    % Extract data from the results structure
    fooof_freqs = fooofresult.freqs; % Frequency values
    power_spectrum = fooofresult.power_spectrum; % Original power spectrum
    % fooofed_spectrum = fooofresult.fooofed_spectrum; % Full FOOOF fit
    ap_fit = fooofresult.ap_fit; % Aperiodic fit
    difference_spectrum = power_spectrum - ap_fit;
    % fooofed_difference_spectrum = fooofed_spectrum - ap_fit;



    % Plot settings
    figure;
    hold on;

    % Plot the original power spectrum
    plot(fooof_freqs, power_spectrum, 'k', 'LineWidth', 1.5, 'DisplayName', 'Power Spectrum');

    % Plot the FOOOFed spectrum (from fooof results)
    % plot(fooof_freqs, fooofed_spectrum, 'r', 'LineWidth', 1.5, 'DisplayName', 'FOOOFed Spectrum');

    % Plot the aperiodic fit
    plot(fooof_freqs, ap_fit, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Aperiodic Fit');

    % Plot the difference spectrum (subtracted in MATLAB)
    plot(fooof_freqs, difference_spectrum, 'g--', 'LineWidth', 1.5, 'DisplayName', 'powerSpect minus apfit');

    % % Plot the FOOOFed difference spectrum (subtracted in MATLAB)
    % plot(fooof_freqs, fooofed_difference_spectrum, '--', 'LineWidth', 1.5, 'DisplayName', 'fooofedSpect minus apfit');

    % Set logarithmic scale for frequency
    % set(gca, 'XScale', 'log'); % Logarithmic x-axis
    % set(gca, 'YScale', 'log'); % Logarithmic y-axis

    % Add labels, legend, and grid
    xlabel('Frequency (Hz)');
    ylabel('Power');
    legend('show');
    grid on;
    title(['S' num2str(SubjG{iStimGroup}(subjectNumInStimGroup)) ' (' FlickerName{iStimGroup} ') - FOOOF Analysis Results: ' num2str(f_range(1)) '-' num2str(f_range(2)) ' Hz']);
    hold off;

end