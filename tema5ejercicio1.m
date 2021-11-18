clear
A = GenerateSparse(7, 7, .4, 52, 115)
B = GenerateSymmetricSparse(5, .8, 52, 115);
[row, col, val] = COO( A );
[rowOff, col, val] = CSR( A );

savingsCOO = zeros(1,100);
savingsCSR = zeros(1,100);

m = 100;
n = 75;
for i=1:101
    density = (i-1)/100;
    % Sean m y n las dimensiones de la matriz original
    % El número de elementos no cero es aproximadamente density*m*n
    % El tamaño de la codificación COO es, por tanto, 3*density*m*n
    % El ahorro sería (3*d*m*n)/(m*n), m*n se cancela y se queda:
    savingsCOO(i) = 3*density;
    % En el caso de CSR, el vector rowOffset tiene tamaño m+1
    savingsCSR(i) = ((m+1)+2*density*m*n)/(m*n);
end

plot(0:0.01:1, savingsCOO);
hold on
plot(0:0.01:1, savingsCSR);
legend('COO', 'CSR');
xlabel('Densidad');
ylabel('Compacta/Original');

function [row, column, value] = COO(matrix)
    compactSize = nnz(matrix);
    row = zeros(1, compactSize);
    column = zeros(1, compactSize);
    value = zeros(1, compactSize);
    index = 1;
    for i=1:size(matrix, 1)
        for j=1:size(matrix, 2)
            if (matrix(i, j) > 1e-30)
                row(index) = i;
                column(index) = j;
                value(index) = matrix(i, j);
                index = index + 1;
            end
        end
    end
end

function [rowOffset, column, value] = CSR(matrix)
    compactSize = nnz(matrix);
    rowOffset = zeros(1, size(matrix, 1) + 1);
    column = zeros(1, compactSize);
    value = zeros(1, compactSize);
    index = 1;
    rowOffset(1)= 0;
    for i=1:size(matrix, 1)
        for j=1:size(matrix, 2)
            if (matrix(i, j) > 1e-30)
                column(index) = j;
                value(index) = matrix(i, j);
                index = index + 1;
            end
        end
        rowOffset(i+1) = index-1;
    end
end

function matrix = GenerateSparse(nRow, nCol, zDensity, vMin, vMax)
    matrix = zeros(nRow, nCol);
    for i=1:nRow
        for j=1:nCol
            if (rand(1) < zDensity)
                matrix(i, j) = vMin + (vMax - vMin) * rand(1);
            end
        end
    end
end

function matrix = GenerateSymmetricSparse(n, zDensity, vMin, vMax)
    matrix = zeros(n);
    for i=1:n
        for j=i:n
            if (rand(1) < zDensity)
                matrix(i, j) = vMin + (vMax - vMin) * rand(1);
            end
        end
    end
    for i=1:n
        for j=1:i
            matrix(i, j) = matrix(j, i);
        end
    end
end