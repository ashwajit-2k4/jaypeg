clear all;       % Clears all variables
clc;             % Clears the command window

% Load a grayscale image
img = imread('barbara256.png');
img = double(img); % Convert to double for DCT calculations
Q = 100;

M = [16 11 10 16 24 40 51 61;
     12 12 14 19 26 58 60 55;
     14 13 16 24 40 57 69 56;
     14 17 22 29 51 87 80 62;
     18 22 37 56 68 109 103 77;
     24 35 55 64 81 104 113 92;
     49 64 78 87 103 121 120 101;
     72 92 95 98 112 100 103 99];

% Define block size (8x8)
blockSize = 8;
[m, n] = size(img);

% Preallocate matrix to store DCT coefficients
dctCoeffs = zeros(m, n);

% Loop through each 8x8 block and compute the DCT
for i = 1:blockSize:m
    for j = 1:blockSize:n
        block = img(i:i+blockSize-1, j:j+blockSize-1);
        dctBlock = dct2(block); % Apply 2D DCT to the block
        dctBlock = round(dctBlock./(M*50/Q));
        dctCoeffs(i:i+blockSize-1, j:j+blockSize-1) = dctBlock;
    end
end

testFrequencies = [1, 2, 3, 4, 5];   % A simple array to check if `sum` works
testSum = sum(testFrequencies);      % This should return 15
disp(testSum);

frequencies = zeros(1,264);
flattenedCoeffs = dctCoeffs(:); % Flatten to a single vector
uniqueValues = unique(flattenedCoeffs);
frequencies = histcounts(flattenedCoeffs, [uniqueValues', max(uniqueValues)+1]);
frequencyTable = table(uniqueValues, frequencies', 'VariableNames', {'Value', 'Frequency'});
% Remove the entry corresponding to 0 from both arrays
size(frequencies)
probabilities = frequencies/sum(frequencies);
symbols = uniqueValues;
dict = huffmandict(symbols, probabilities);
disp(dict)

% Example 8x8 matrix
A = dctCoeffs(249:256,249:256); % This creates a matrix with values from 1 to 64

% Initialize the result array
zigzagArray = zeros(64);
i = 1;
% Zigzag traversal
index = 1;
flag = 1;
for sum = 2:16
    if flag == 0
        for row = 1:8
            if row < sum && sum-row <= 8
                zigzagArray(i) = A(row,sum-row);
                disp(zigzagArray(i));
                i = i+1;
            end
        end
        flag = 1;
    else
        for row = 8:-1:1
            if row < sum && sum-row <= 8
                zigzagArray(i) = A(row,sum-row);
                disp(zigzagArray(i));
                i = i+1;
            end
        end
        flag = 0;
    end
end

% Display the result
disp(dctCoeffs(249:256,249:256))


