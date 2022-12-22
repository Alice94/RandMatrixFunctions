figure('Position', [100 100 600 300]);
Legend = [];

load(strcat(NAME, ".mat"));

if (algorithms(1) == 1)
    semilogy(ms, min(1, errArnoldi), '-ko', 'linewidth', 1)
    hold on
    Legend = [Legend; 'Arnoldi '];
end

if (algorithms(2) == 1)
    semilogy(ms, min(1, errBG), '-m*', 'linewidth', 1)
    hold on
    Legend = [Legend; 'BG-f(A) '];
end

if (algorithms(3) == 1)
    semilogy(ms, min(1, errNTB), '-c.', 'linewidth', 1)
    hold on
    Legend = [Legend; 'NT1-f(A)'];
end

if (algorithms(4) == 1)
    semilogy(ms, min(1, errNTtrid), '-x', 'linewidth', 1, 'color', "#EDB120")
    hold on
    Legend = [Legend; 'NT2-f(A)'];
end

if (algorithms(5) == 1)
    semilogy(ms, min(1, errsFOM), '-rh', 'linewidth', 1)
    hold on
    Legend = [Legend; 'sFOM    '];
end

if (algorithms(6) == 1)
    semilogy(ms, min(1, errNTW), '--b^', 'linewidth', 1)
    hold on
    Legend = [Legend; "Whiten  "];
end



xlabel('m');
ylabel('relative error');
legend(Legend, 'Location', 'best')
hold off

saveas(gcf,NAME,'epsc')