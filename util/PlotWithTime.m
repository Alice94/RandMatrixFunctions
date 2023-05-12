figure('Position', [100 100 800 300]);
Legend1 = [];

load(strcat(NAME, ".mat"));

%%

if (algorithms(1) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errArnoldi), '-ko', 'linewidth', 1)
    hold on
    subplot(1,2,2)
    plot(ms, timeArnoldi, '-ko', 'linewidth', 1)
    hold on
    Legend1 = [Legend1; 'Arnoldi '];
end

if (algorithms(2) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errBG), '-m*', 'linewidth', 1)
    hold on
    subplot(1,2,2)
    plot(ms, timeBG, '-m*', 'linewidth', 1)
    hold on
    Legend1 = [Legend1; 'BG-f(A) '];
end

if (algorithms(3) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errNTB), '-c.', 'linewidth', 1)
    hold on
    subplot(1,2,2)
    plot(ms, timeNTB, '-c.', 'linewidth', 1)
    hold on
    Legend1 = [Legend1; 'NT1-f(A)'];
end

if (algorithms(4) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errNTtrid), '-x', 'linewidth', 1, 'color', "#EDB120")
    hold on
    subplot(1,2,2)
    plot(ms, timeNTtrid, '-x', 'linewidth', 1, 'color', "#EDB120")
    hold on
    Legend1 = [Legend1; 'NT2-f(A)'];
end

if (algorithms(6) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errNTW), '--b^', 'linewidth', 1)
    hold on
    subplot(1,2,2)
    plot(ms, timeNTW, '--b^', 'linewidth', 1)
    hold on
    Legend1 = [Legend1; "Whiten  "];
end

if (algorithms(7) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errArnoldiRestart), '-g*', 'linewidth', 1)
    hold on
    subplot(1,2,2)
    plot(ms, timeArnoldiRestart, '-g*', 'linewidth', 1)
    hold on
    Legend1 = [Legend1; "Restart "];
end

Legend2 = Legend1;

if (algorithms(5) == 1)
    subplot(1,2,1)
    semilogy(ms, min(1, errsFOM), '-rh', 'linewidth', 1)
    hold on
    subplot(1,2,2)
    plot(ms, timesFOM, '-rh', 'linewidth', 1)
    hold on 
    Legend1 = [Legend1; 'sFOM    '];
end

subplot(1,2,1)
xlabel('m');
ylabel('relative error');
subplot(1,2,2)
xlabel('m');
ylabel('time (seconds)')

subplot(1,2,1)
legend(Legend1, 'Location', 'best')
hold off
subplot(1,2,2)
legend(Legend1, 'Location', 'best')
hold off

saveas(gcf,NAME,'epsc')