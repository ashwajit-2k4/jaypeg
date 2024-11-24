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

% Flatten the DCT coefficients to a single vector
flattenedCoeffs = dctCoeffs(:);

% Display the modified flattenedCoeffs
%disp(flattenedCoeffs);

% Get unique values and their frequencies
uniqueValues = unique(flattenedCoeffs);
frequencies = histcounts(flattenedCoeffs, [uniqueValues', max(uniqueValues)+1]);

% Remove the entry corresponding to 0 from both arrays
nonZeroMask = uniqueValues ~= 0; % Mask to exclude 0
uniqueValues = uniqueValues(nonZeroMask);  % Exclude 0 from unique values
frequencies = frequencies(nonZeroMask);   % Exclude 0's frequency

% Create the frequency table
frequencyTable = table(uniqueValues, frequencies', 'VariableNames', {'Value', 'Frequency'});

% Calculate probabilities (excluding 0)
probabilities = frequencies / sum(frequencies);

% Define the symbols (without 0)
symbols = uniqueValues;

% Create the Huffman dictionary
dict = huffmandict(symbols, probabilities);

% Display the Huffman dictionary
% disp("Huffman Dictionary:");
% disp(dict);

% Example 8x8 matrix
A = dctCoeffs(249:256,249:256); % This creates a matrix with values from 1 to 64

% Initialize the result array
zigzagArray = zeros(64);
i = 1;
% Zigzag traversal
index = 1;
flag = 1;
width = 256;
height = 256;
dcCoeffs = zeros(1,1024);
dcCoeffs2 = zeros(1,1024);
prev_zig = 0;

for k = 1:(width/8)
    for j = 1:(height/8)
        flag = 1;
        i=1;
        minx=(k-1)*8+1;
        maxx = 8*k;
        miny = (j-1)*8+1;
        maxy = 8*j;
        A = dctCoeffs(minx:maxx,miny:maxy);
        zigzagArray = zeros(64,1);
        for sum = 2:16
            if flag == 0
                for row = 1:8
                    if row < sum && sum-row <= 8
                        zigzagArray(i) = A(row,sum-row);
                        i = i+1;
                    end
                end
                flag = 1;
            else
                for row = 8:-1:1
                    if row < sum && sum-row <= 8
                        zigzagArray(i) = A(row,sum-row);
                        i = i+1;
                     
                    end
                end
                flag = 0;
            end
        end
        orderedCoeffs(1,2016*(j-1)+63*(k-1)+1:2016*(j-1)+63*k) = zigzagArray(2:64,1);
        dcCoeffs(32*(j-1)+k) = zigzagArray(1,1) - prev_zig;
        dcCoeffs2(32*(j-1)+k) = zigzagArray(1,1);

        prev_zig = zigzagArray(1,1);
    end
end
% disp(orderedCoeffs)
zigzagArray = zigzagArray(2:end);
% Display the result
% disp(dctCoeffs(249:256,249:256))

imageRun = orderedCoeffs;

results = []; % Each row will contain: [Non-zero Number, Zero Count, Code Length, Huffman Code]
% Initialize index for traversing the run
i = 1;
results = []; % Initialize results
i = 1; % Starting index

while i <= length(imageRun)
    if imageRun(i) ~= 0
        % Non-zero number found
        nonZeroNum = imageRun(i);
        
        % Count consecutive zeroes
        zeroCount = 0;
        j = i + 1;
        while j <= length(imageRun) && imageRun(j) == 0
            zeroCount = zeroCount + 1;
            j = j + 1;
        end
        
        % Calculate Huffman code length (ceil(log2(abs(nonZeroNum) + 1)))
        codeLength = ceil(log2(abs(nonZeroNum) + 1));
        
        % Generate Huffman code (binary representation)
        huffmanCode = dec2bin(nonZeroNum, codeLength); % Binary representation
        index = find(symbols == imageRun(i));
        
        % Retrieve the corresponding Huffman code for the symbol
        code = dict{index, 2};
        huffmanCode = num2str(code); % Ensure it's a string
        huffmanCode = strrep(huffmanCode, ' ', ''); % Remove any spaces (optional)
        
        % Store the result for the non-zero number
        results = [results; {nonZeroNum, min(zeroCount - floor(zeroCount / 15) * 15, 15), length(huffmanCode), huffmanCode}];
        
        % Handle large zero counts (> 15)
        zeroCount = zeroCount - 15;
        while zeroCount >= 0
            % Add an entry for a chunk of 15 zeroes
            results = [results; {0, 15, 0, '0'}];
            zeroCount = zeroCount - 15;
        end
        
        % Add the remaining zero count if > 0
        if zeroCount > 0
            results = [results; {0, zeroCount, 0, '0'}];
        end
        
        % Move index forward
        i = j; % Skip processed elements
    else
        i = i + 1; % Skip zero
    end
end

% Display results
%disp('Results:');
disp('Non-Zero | Zero Count | Code Length | Huffman Code');
disp(results);


% Assuming results is a 2D array with columns:
% Column 1: Non-zero number
% Column 2: Zero count
% Column 3: Code length
% Column 4: Huffman code (binary string)

binaryArray = {};  % Cell array to hold the binary strings for each entry

% Loop through each row of the results array
binaryString = '';
for i = 1:size(results)
    % 1. Get the zero count (4 bits)
    zeroCount = results{i, 2};
    zeroCountBits = dec2bin(zeroCount, 4); % Convert to 4-bit binary string
    %disp(zeroCountBits)
    % 2. Get the size of the Huffman code (code length) (4 bits)
    codeLength = results{i, 3};
    codeLengthBits = dec2bin(codeLength, 4); % Convert to 4-bit binary string
    %disp(codeLengthBits)
    % 3. Get the Huffman code (binary string)
    huffmanCode = results{i, 4}; % Extract the Huffman code from the table
    disp(huffmanCode)
    % 4. Concatenate all parts: zero count (4 bits), code length (4 bits), and Huffman code
    binaryString = [binaryString, zeroCountBits, codeLengthBits, huffmanCode];
    
    % 5. Store the binary string in the array
    binaryArray = binaryString;
end
encodedAc = binaryString;
% Display the resulting binary array
% disp('Binary Array:');
% disp(binaryArray)
% disp(length(binaryArray))
% disp(dcCoeffs)
encodedDc = '';
for i = 1:length(dcCoeffs2)
    index = find(symbols == dcCoeffs2(i));
    %disp(index)
% Retrieve the corresponding Huffman code for the symbol
    code = dict{index, 2};
    %disp(code)
    binaryString = dec2bin(code);
    %disp(binaryString')
    encodedDc = [encodedDc, binaryString'];
    % Remove any spaces between the digits (optional)
end
%disp(encodedDc)

results2 = {};

% Loop through each symbol and code
encodedHuff = '';
%disp("hi")
for idx = 1:length(symbols)
    num = dict{idx,1};
    code = dict{idx,2};
    % a) Calculate the number of bits required to represent the number (signed, sign magnitude)
    numBits = floor(log2(abs(num))) + 1 + 1;  % Extra bit for sign magnitude (sign bit)
    
    % b) Store the number of bits required as an 8-bit binary number
    numBitsBinary = dec2bin(numBits, 4);  % Store as 8-bit binary
    
    % c) The number itself in sign magnitude binary representation (with sign bit at the bottom)
    if num < 0
        % Convert negative number to sign magnitude representation
        magnitudeBinary = dec2bin(abs(num), numBits - 1);  % Get the magnitude in binary (without sign bit)
        %disp(magnitudeBinary)
        numberBinary = [magnitudeBinary, '1'];  % Append the sign bit (1) at the bottom (LSB)
    else
        numberBinary = [dec2bin(num, numBits - 1), '0'];  % Positive number, append '0' as sign bit at the bottom
    end
    
    % d) Number of bits in the Huffman code
    huffmanBits = length(code);
    huffmanBitsBinary = dec2bin(huffmanBits, 4);  % Convert to binary (8-bit format)
    
    % e) The Huffman code itself
    huffmanCode = num2str(code);
    huffmanCode = strrep(huffmanCode, ' ', ''); % Remove spaces (e.g., '10100110010110')
    %disp(huffmanCode)

    % Store the result in binary format
    results2{idx} = sprintf('%s   %s   %s   %s', ...
    numBitsBinary, numberBinary, huffmanBitsBinary, huffmanCode);
    encodedHuff = [encodedHuff, numBitsBinary, numberBinary, huffmanBitsBinary, huffmanCode];
end

% Binary strings to write
binaryStrings = [encodedHuff, encodedDc, encodedAc];
% Binary strings as rows

% Open a file for writing
fileID = fopen('example2.txt', 'w'); % 'w' for write mode (overwrites if the file exists)

% Data to write
data = [1, 2, 3; 4, 5, 6];

% Write each row of the array to the file
fprintf(fileID, '%s\n', encodedHuff);
fprintf(fileID, '%s\n', encodedDc);
fprintf(fileID, '%s\n', encodedAc);


% Close the file
fclose(fileID);
