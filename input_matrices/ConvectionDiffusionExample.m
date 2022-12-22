function A = ConvectionDiffusionExample(N, Pe)

h = 1/(N-1);
T = spdiags(ones(N-2, 1) * [-1, 2, -1], -1:1, N-2, N-2);
A1 = kron(speye(N-2), T) * 1/2;
A2 = kron(T, speye(N-2));

tic;
for i = 1:N-2
    for j = 1:N-2
        if (i/(N-1) >= 0.25 && i/(N-1) <= 0.75 && j/(N-1) >= 0.25 && j/(N-1) <= 0.75)
            A1((i-1)*(N-2) + j, (i-1)*(N-2) + j) = 1000*A1((i-1)*(N-2) + j, (i-1)*(N-2) + j);
            if ((j+1)/(N-1) <= 0.75)
                A1((i-1)*(N-2) + j, (i-1)*(N-2) + j+1) = ...
                    1000*A1((i-1)*(N-2) + j, (i-1)*(N-2) + j+1);
            end
            if ((j-1)/(N-1) >= 0.25)
                A1((i-1)*(N-2) + j, (i-1)*(N-2) + j-1) = ...
                    1000*A1((i-1)*(N-2) + j, (i-1)*(N-2) + j-1);
            end
            
            A2((i-1)*(N-2) + j, (i-1)*(N-2) + j) = 1000*A2((i-1)*(N-2) + j, (i-1)*(N-2) + j);
            if ((i+1)/(N-1) <= 0.75)
                A2((i-1)*(N-2) + j, i*(N-2) + j) = ...
                    1000*A2((i-1)*(N-2) + j, i*(N-2) + j);
            end
            if ((i-1)/(N-1) >= 0.25)
                A2((i-1)*(N-2) + j, (i-2)*(N-2) + j) = ...
                    1000*A2((i-1)*(N-2) + j, (i-2)*(N-2) + j);
            end
        end
    end
end
toc
tic;
A3 = sparse((N-2)^2, (N-2)^2);
for i = 1:N-2
    for j = 1:N-2
        if i < N-2, A3((i-1)*(N-2) + j, i*(N-2) + j) = 2*i + 2*j + 1; end
        if j < N-2, A3((i-1)*(N-2) + j, (i-1)*(N-2) + j + 1) = 2*i - 2*j - 1; end
    end
end
A3 = A3 - A3';
A3 = A3 * Pe / 4;
toc

A = A1 + A2 + A3*h^2;

disp(normest((A+A')/2, 1e-2))
disp(normest((A-A')/2, 1e-2))