clc; clear; close all;

choice = menu('Select Animation', ...
    '1) rect * rect -> triangle (improved with labels and colors)', ...
    '2) Periodic triangular message (improved with cumulative build-up)', ...
    '3) Building the spectrum of periodic triangular message', ...
    '4) AM modulation envelope (improved true animation with zoom and labels)', ...
    '5) Effect of modulation index in AM (from under to over-modulation)', ...
    '6) AM modulated spectrum animation (shifting message spectrum)', ...
    '7) DSB-SC phase effect (improved with simultaneous animation and envelopes)', ...
    '8) Continuous phase change in DSB-SC', ...
    '9) Comparison between AM and DSB-SC waveforms');

switch choice
%% =========================================================
case 1 % rect * rect -> triangle (improved with labels, colors, and slower animation)
    dt = 0.002;
    t = -2:dt:2;
    rect1 = double(abs(t) <= 0.5);
    figure('Color','w', 'Position', [100 100 800 600])
    subplot(2,1,1)
    h1 = plot(t, rect1, 'b', 'LineWidth', 2); hold on
    h2 = plot(t, zeros(size(t)), 'r', 'LineWidth', 2);
    h3 = area(t, zeros(size(t)), 'FaceColor', [0.8 0.8 1], 'EdgeColor','none');
    title('Convolution: Sliding Rectangles', 'FontSize', 14)
    xlabel('Time t', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-0.1 1.2]); grid on
    legend('Fixed rect(t)', 'Sliding rect(t - \tau)', 'Overlap Area')
    subplot(2,1,2)
    h4 = animatedline('Color', 'g', 'LineWidth', 2);
    title('Resulting Triangle from Convolution', 'FontSize', 14)
    xlabel('Shift \tau', 'FontSize', 12)
    ylabel('Convolution Value', 'FontSize', 12)
    ylim([0 1.2]); grid on
    for k = 1:length(t)
        shift = t(k);
        rect2 = double(abs(t - shift) <= 0.5);
        overlap = rect1 .* rect2;
        conv_val = sum(overlap)*dt;
        set(h2,'YData',rect2)
        set(h3,'YData',overlap)
        addpoints(h4, t(k), conv_val)
        drawnow
        pause(0.005) % Slower for better visualization
    end
%% =========================================================
case 2 % Periodic triangular message (improved with cumulative build-up of pulses)
    dt = 0.0005;
    t = -0.4:dt:0.4;
    T0 = 0.12;
    alpha = 0.04;
    figure('Color','w', 'Position', [100 100 800 400])
    h = animatedline('LineWidth', 2, 'Color', 'b');
    title('Building Periodic Triangular Message m(t)', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([0 1.2]); grid on
    m_cum = zeros(size(t));
    for n = -5:5
        tau = t - n*T0;
        tri = max(1 - abs(tau)/alpha, 0);
        m_cum = m_cum + tri;
        set(h, 'XData', t, 'YData', m_cum);
        drawnow
        pause(0.5) % Pause to show each pulse addition
    end
    text(0, 1.1, 'Final Periodic Signal', 'FontSize', 12, 'HorizontalAlignment', 'center')
%% =========================================================
case 3 % Building the spectrum of periodic triangular message
    alpha = 40e-3;
    f0 = 1 / 120e-3;
    n_max = 10; % Number of harmonics to show
    f = -n_max*f0 : f0/100 : n_max*f0;
    sinc_env = alpha * f0 * sinc(alpha * f).^2; % Fixed: added f0
    figure('Color','w', 'Position', [100 100 800 400])
    plot(f, sinc_env, 'k--', 'LineWidth', 1.5); hold on
    title('Spectrum of Periodic Triangular m(t): Building Harmonics', 'FontSize', 14)
    xlabel('Frequency (Hz)', 'FontSize', 12)
    ylabel('|M(f)|', 'FontSize', 12)
    grid on
    ylim([0 max(sinc_env)*1.1])
    h_stems = [];
    for n = -n_max:n_max
        cn = alpha * f0 * sinc(n * alpha * f0)^2;
        h_stem = stem(n*f0, cn, 'b', 'LineWidth', 2, 'MarkerSize', 8);
        h_stems = [h_stems h_stem];
        text(n*f0, cn + 0.01, ['c_' num2str(n)], 'FontSize', 10, 'HorizontalAlignment', 'center')
        drawnow
        pause(0.5)
    end
    legend('sinc^2 Envelope', 'Harmonic Components')
%% =========================================================
case 4 % AM modulation TRUE animation (improved with moving window, labels, and envelope highlights)
    dt = 0.0005;
    t = 0:dt:0.4;
    fc = 100;
    T0 = 0.12;
    alpha = 0.04;
    % build periodic triangular message
    m = zeros(size(t));
    for n = -5:5
        tau = t - n*T0;
        tri = max(1 - abs(tau)/alpha, 0);
        m = m + tri;
    end
    carrier = cos(2*pi*fc*t);
    s = (1 + m) .* carrier; % Assuming ka=1, Ac=1
    window = 0.03; % visible moving time window
    figure('Color','w', 'Position', [100 100 800 800])
    subplot(3,1,1)
    plot(t, m, 'b', 'LineWidth', 2)
    title('Message Signal m(t)', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([0 1.2]); grid on
    h_cursor1 = xline(0,'r','LineWidth',2);
    subplot(3,1,2)
    plot(t, carrier, 'r', 'LineWidth', 2)
    title('Carrier Signal c(t)', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-1.2 1.2]); grid on
    h_cursor2 = xline(0,'r','LineWidth',2);
    subplot(3,1,3)
    plot(t, s, 'g', 'LineWidth', 2); hold on
    h_env_pos = plot(t, 1+m, 'k--', 'LineWidth', 1.5);
    h_env_neg = plot(t, -(1+m), 'k--', 'LineWidth', 1.5);
    title('AM Modulated Signal s(t) with Envelopes', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-2 2]); grid on
    legend('s(t)', 'Positive Envelope', 'Negative Envelope')
    h_cursor3 = xline(0,'r','LineWidth',2);
    for k = 1:5:length(t)
        t_now = t(k);
        set(h_cursor1,'Value',t_now)
        set(h_cursor2,'Value',t_now)
        set(h_cursor3,'Value',t_now)
        subplot(3,1,1)
        xlim([t_now-window/2 t_now+window/2]) % Centered zoom
        subplot(3,1,2)
        xlim([t_now-window/2 t_now+window/2])
        subplot(3,1,3)
        xlim([t_now-window/2 t_now+window/2])
        drawnow
        pause(0.01)
    end
%% =========================================================
case 5 % Effect of modulation index in AM (improved with triangular message)
    dt = 0.0005;
    t = 0:dt:0.4; % Longer to show multiple periods
    fc = 100; % Higher fc for better visualization
    T0 = 0.12;
    alpha = 0.04;
    Ac = 1;
    % build periodic triangular message, max=1
    m = zeros(size(t));
    for n = -5:5
        tau = t - n*T0;
        tri = max(1 - abs(tau)/alpha, 0);
        m = m + tri;
    end
    figure('Color','w', 'Position', [100 100 800 400])
    h_s = animatedline('Color', 'b', 'LineWidth', 2);
    hold on
    h_env_pos = plot(t, zeros(size(t)), 'k--', 'LineWidth', 1.5);
    h_env_neg = plot(t, zeros(size(t)), 'k--', 'LineWidth', 1.5);
    title('AM with Varying Modulation Index \mu (Triangular Message)', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-2 2]); grid on
    h_text = text(0.05, 1.5, '\mu = 0', 'FontSize', 12);
    mu_values = 0:0.05:1.5; % Finer steps
    for mu = mu_values
        s = Ac * (1 + mu * m) .* cos(2*pi*fc*t);
        env_pos = Ac * (1 + mu * m);
        env_neg = -Ac * (1 + mu * m);
        set(h_env_pos, 'YData', env_pos);
        set(h_env_neg, 'YData', env_neg);
        set(h_text, 'String', ['\mu = ' num2str(mu)]);
        set(h_s, 'XData', [], 'YData', []);
        for k = 1:length(t)
            xdata = get(h_s, 'XData');
            ydata = get(h_s, 'YData');
            set(h_s, 'XData', [xdata t(k)], 'YData', [ydata s(k)]);
            drawnow limitrate
        end
        pause(0.5)
    end
    text(0.05, -1.5, 'Over-modulation causes envelope distortion', 'FontSize', 10)
%% =========================================================
case 6 % AM modulated spectrum animation (improved with discrete harmonics)
    alpha = 40e-3;
    f0 = 1 / 120e-3;
    fc = 100; % Carrier frequency
    n_max = 5;
    n = -n_max:n_max;
    f_disc = n * f0;
    M_disc = alpha * f0 * sinc(alpha * f_disc).^2;
    figure('Color','w', 'Position', [100 100 800 400])
    subplot(2,1,1)
    stem(f_disc, abs(M_disc), 'b', 'LineWidth', 2); hold on
    title('Message Spectrum M(f) - Discrete Harmonics', 'FontSize', 14)
    xlabel('Frequency (Hz)', 'FontSize', 12)
    ylabel('|M(f)|', 'FontSize', 12)
    grid on
    subplot(2,1,2)
    title('AM Spectrum S(f): Shifting and Adding', 'FontSize', 14)
    xlabel('Frequency (Hz)', 'FontSize', 12)
    ylabel('|S(f)|', 'FontSize', 12)
    xlim([-fc - n_max*f0 - 10, fc + n_max*f0 + 10])
    ylim([0 max([0.5, max(0.5*abs(M_disc))])*1.1])
    grid on
    % Animate shifting to +fc and -fc
    for shift = 0:2:fc % Larger step for faster animation
        cla;
        stem([-fc fc], [0.5 0.5], 'r', 'LineWidth', 2);
        f_pos = f_disc + shift;
        f_neg = f_disc - shift;
        stem(f_pos, 0.5 * abs(M_disc), 'g');
        stem(f_neg, 0.5 * abs(M_disc), 'g');
        xlim([-fc - n_max*f0 - 10, fc + n_max*f0 + 10])
        ylim([0 max([0.5, max(0.5*abs(M_disc))])*1.1])
        grid on
        drawnow
        pause(0.01)
    end
    legend('Carrier Deltas', 'Shifted M(f)')
%% =========================================================
case 7 % DSB-SC phase effect (improved with simultaneous point-by-point animation across all phases)
    t = 0:0.0005:0.2;
    fm = 5;
    fc = 50;
    phi_list = [0 pi/4 pi/2 3*pi/4];
    Am = 1; Ac = 1;
    figure('Color','w', 'Position', [100 100 800 800])
    h = gobjects(length(phi_list),1);
    s_list = zeros(length(phi_list), length(t));
    m = Am * cos(2*pi*fm*t);
    for p = 1:length(phi_list)
        phi = phi_list(p);
        subplot(4,1,p)
        h(p) = animatedline('Color', 'm', 'LineWidth', 2);
        hold on
        plot(t, m, 'b--', 'LineWidth', 1) % Message envelope
        plot(t, -m, 'b--', 'LineWidth', 1)
        ylim([-1.2 1.2]); grid on
        title(['DSB-SC with Phase \phi = ' num2str(phi*180/pi) '°'], 'FontSize', 12)
        xlabel('Time (s)', 'FontSize', 10)
        ylabel('Amplitude', 'FontSize', 10)
        legend('s(t)', 'Envelope', '')
        c = Ac * cos(2*pi*fc*t + phi);
        s_list(p,:) = m .* c;
    end
    for k = 1:length(t)
        for p = 1:length(phi_list)
            xdata = get(h(p), 'XData');
            ydata = get(h(p), 'YData');
            set(h(p), 'XData', [xdata t(k)], 'YData', [ydata s_list(p,k)]);
        end
        drawnow limitrate
    end
%% =========================================================
case 8 % Continuous phase change in DSB-SC
    t = 0:0.0005:0.2;
    fm = 5;
    fc = 50;
    Am = 1; Ac = 1;
    figure('Color','w', 'Position', [100 100 800 400])
    h_s = animatedline('Color', 'c', 'LineWidth', 2);
    title('DSB-SC with Continuous Phase Change', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-1.2 1.2]); grid on
    h_text = text(0.1, 1.0, '\phi = 0°', 'FontSize', 12);
    m = Am * cos(2*pi*fm*t);
    for phi_deg = 0:5:180
        phi = phi_deg * pi / 180;
        c = Ac * cos(2*pi*fc*t + phi);
        s = m .* c;
        set(h_s, 'XData', [], 'YData', []);
        for k = 1:length(t)
            xdata = get(h_s, 'XData');
            ydata = get(h_s, 'YData');
            set(h_s, 'XData', [xdata t(k)], 'YData', [ydata s(k)]);
            drawnow limitrate
        end
        set(h_text, 'String', ['\phi = ' num2str(phi_deg) '°'])
        pause(0.5)
    end
%% =========================================================
case 9 % Comparison between AM and DSB-SC waveforms (fixed to animate without initial full plot)
    dt = 0.0005;
    t = 0:dt:0.2;
    fm = 5;
    fc = 50;
    Am = 1; Ac = 1;
    phi = 0; % For DSB-SC
    m = Am * cos(2*pi*fm*t);
    c = Ac * cos(2*pi*fc*t + phi);
    s_dsb = m .* c;
    s_am = (1 + m) .* c; % mu=1
    figure('Color','w', 'Position', [100 100 800 600])
    subplot(2,1,1)
    plot(t, 1 + m, 'k--'); hold on
    plot(t, -(1 + m), 'k--')
    title('AM Waveform (with Carrier)', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-2 2]); grid on
    legend('Envelopes')
    h_am = animatedline('Color','b', 'LineWidth', 2);
    subplot(2,1,2)
    plot(t, m, 'k--'); hold on
    plot(t, -m, 'k--')
    title('DSB-SC Waveform (no Carrier)', 'FontSize', 14)
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    ylim([-1.2 1.2]); grid on
    legend('Envelopes')
    h_dsb = animatedline('Color','r', 'LineWidth', 2);
    for k = 1:length(t)
        xdata_am = get(h_am, 'XData');
        ydata_am = get(h_am, 'YData');
        set(h_am, 'XData', [xdata_am t(k)], 'YData', [ydata_am s_am(k)]);
        xdata_dsb = get(h_dsb, 'XData');
        ydata_dsb = get(h_dsb, 'YData');
        set(h_dsb, 'XData', [xdata_dsb t(k)], 'YData', [ydata_dsb s_dsb(k)]);
        drawnow limitrate
    end
end