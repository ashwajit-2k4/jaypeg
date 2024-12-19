clc;
quality = 10:15:100;
num_quality_factors = length(quality);
image_files = dir('bw/*.png');
num_images = length(image_files);

RMSE = zeros(num_images, num_quality_factors);
BPP = zeros(num_images, num_quality_factors);

for a = 2:2
    clearvars -except a image_files num_images quality num_quality_factors; % Clears all variables
    filename= sprintf('bw/%d.png', a);
    img = imread(filename);
    [height, width] = size(img);
    num_pixels = height * width;

    for b = 1:7
        Q = quality(b);

       % encoded_image = encode_jpeg(original_image, quality);
       % img = imread('kodak242.png');
img = double(img); % Convert to double for DCT calculations

% The quantisation matrix and the corresponding scaling factor

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
%[height, width] = size(img);

% The amount of mean padding to make height and width a multiple of 8
padHeight = 8-mod(height, 8);
padWidth = 8-mod(width, 8);

% If padding is needed, calculate the mean and pad
if padHeight ~= 8
    meanRow = mean(img(end-(8-padHeight)+1:end, :), 1); % Mean of the last row
    img = [img; repmat(meanRow, padHeight, 1)]; % Add rows with the mean value
end

if padWidth ~= 8
    meanCol = mean(img(:, end-(8-padWidth)+1:end), 2); % Mean of the last columns
    img = [img, repmat(meanCol, 1, padWidth)]; % Add columns with the mean value
end

% Verify new dimensions are multiples of 8
[height, width] = size(img);
disp(['New height: ', num2str(height)]);
disp(['New width: ', num2str(width)]);  
dctCoeffs = []
% Stores the newly calculated DCT coefficients
dctCoeffs = zeros(height, width);

% Loop through each 8x8 block and compute the DCT
for i = 1:blockSize:(height-7)
    for j = 1:blockSize:(width-7)
        
        % Selects an 8x8 block
        block = img(i:i+blockSize-1, j:j+blockSize-1);
        dctBlock = dct2(block); % Apply 2D DCT to the block
        dctBlock = round(dctBlock./(M*50/Q)); % Quantise by element wise division by M/Q
        dctCoeffs(i:i+blockSize-1, j:j+blockSize-1) = dctBlock; % Assign back to the DCT matrix
    end
end

numBlocksHeight = height / 8;   % Number of 8x8 blocks in the height direction
numBlocksWidth = width / 8;     % Number of 8x8 blocks in the width direction

% Initialize an array to store the encoded DC elements
encodedTopLeft = zeros(numBlocksHeight, numBlocksWidth);

% Initialize the previous top-left element (for differential encoding)
% since the first element is stored as is
previousTopLeft = 0; 

% The zero count of DC coefficient differences (since they need a Huffman
% code, but the AC ones do not)
zcount = 0;

% Loop over all blocks (top-left to bottom-right)
for blockCol = 1:numBlocksWidth
    for blockRow = 1:numBlocksHeight
        % Define the block's row and column indices
        minRow = (blockRow - 1) * 8 + 1;
        minCol = (blockCol - 1) * 8 + 1;
        
        % Extract the DC (top-left) element of the current block
        currentTopLeft = dctCoeffs(minRow, minCol);
        
        % Store elements as its difference with the previous
        encodedTopLeft(blockRow, blockCol) = currentTopLeft - previousTopLeft;
        
        % Update the previous top-left element
        previousTopLeft = currentTopLeft;

        % Assign it back to dctcoeffs
        dctCoeffs(minRow, minCol) = encodedTopLeft(blockRow,blockCol);

        % Check for zero count
        if encodedTopLeft(blockRow,blockCol) == 0
            zcount = zcount+1;
        end
    end
end

% Flatten the DCT coefficients to a single vector
flattenedCoeffs = dctCoeffs(:);
% Get unique values and their frequencies for Huffman encoding

uniqueValues = unique(flattenedCoeffs);


frequencies = histcounts(flattenedCoeffs, [uniqueValues', max(uniqueValues)+1]);

% Create the frequency table
frequencyTable = table(uniqueValues, frequencies', 'VariableNames', {'Value', 'Frequency'});
% Finds the index corresponding to 0
zeroloc = find(uniqueValues == 0);
% Only stores the DC coefficients 0
frequencies(zeroloc) = zcount; 
% Calculate probabilities of each number occurring
x = 0;
for i = 1:length(frequencies)
    x = x+frequencies(i);
end
probabilities = frequencies / x;
% Create the Huffman dictionary
dict = huffmandict(uniqueValues, probabilities);

% Display the Huffman dictionary
disp("Huffman Dictionary:");
disp(dict);

% For arranging elements in zigzag order
% Initialize the result array
zigzagArray = zeros(64);

% Used to switch directions at each end
flag = 1;
dcCoeffs = zeros(1,height*width/64);
% Iterate over all blocks
for k = 1:(width/8)
    for j = 1:(height/8)
        
        % Set flag, i, start and end points
        flag = 1;
        i=1;
        minx=(k-1)*8+1;
        maxx = 8*k;
        miny = (j-1)*8+1;
        maxy = 8*j;

        % Take a DCT block
        A = dctCoeffs(miny:maxy,minx:maxx);
        zigzagArray = zeros(64,1);

        % Follow zigzag pattern, flag allows change of direction
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

        % Stores AC coefficients
        acCoeffs(1,height*8*63/64*(k-1)+63*(j-1)+1:height*8*63/64*(k-1)+63*j) = zigzagArray(2:64,1);

        % Stores DC coefficients
        dcCoeffs(height/8*(k-1)+j) = zigzagArray(1,1);
    end
end

results = []; % Initialize results
i = 1; % Starting index

% Iterates over AC coefficients
while i <= length(acCoeffs)
    % Checks if non-zero
    if acCoeffs(i) ~= 0

        nonZeroNum = acCoeffs(i);
        
        % Count consecutive zeroes
        zeroCount = 0;
        j = i + 1;
        while j <= length(acCoeffs) && acCoeffs(j) == 0
            zeroCount = zeroCount + 1;
            j = j + 1;
        end
                
        % Find index corresponding to the number
        index = find(uniqueValues == acCoeffs(i));
        % Retrieve the corresponding Huffman code for the symbol
        code = dict{index, 2};
        huffmanCode = num2str(code); % Ensure it's a string
        huffmanCode = strrep(huffmanCode, ' ', ''); % Remove any spaces
        
        % Store the result for the non-zero number, stores zeroCount - 15n,
        % rest stored as (15,0,0)
        results = [results; {nonZeroNum, min(zeroCount - floor(zeroCount / 15) * 15, 15), length(huffmanCode), huffmanCode}];
        
        % Handle large zero counts (> 15)
        zeroCount = zeroCount - 15;
        while zeroCount >= 0
            % Add an entry for a chunk of 15 zeroes
            results = [results; {0, 15, 0, '0'}];
            zeroCount = zeroCount - 15;
        end
        
        % Adds the remaining zero count if it is there. Not strictly
        % necessary since we account for it previously
        if zeroCount > 0
            results = [results; {0, zeroCount, 0, '0'}];
        end
        
        i = j; % Skip processed elements
    else
        i = i + 1; % Skip zero
    end
end

disp('Non-Zero | Zero Count | Code Length | Huffman Code');
disp(results);

% Loop through each row of the results array
encodedAc = '';
for i = 1:size(results)

    % Get the zero count (4 bits)
    zeroCount = results{i, 2};
    zeroCountBits = dec2bin(zeroCount, 4); % Convert to 4-bit binary string

    % 2. Get the size of the Huffman code (code length) (6 bits)
    codeLength = results{i, 3};
    codeLengthBits = dec2bin(codeLength, 6); % Convert to 6-bit binary string
    
    % 3. Get the Huffman code (binary string)
    huffmanCode = results{i, 4}; % Extract the Huffman code from the table

    % 4. Concatenate all parts: zero count (4 bits), code length (6 bits), and Huffman code
    encodedAc = [encodedAc, zeroCountBits, codeLengthBits, huffmanCode];
end

encodedDc = '';
for i = 1:length(dcCoeffs)

    % Finds index corresponding to number
    index = find(uniqueValues == dcCoeffs(i));

    % Retrieve the corresponding Huffman code for the symbol
    code = dict{index, 2};
    binaryString = dec2bin(code); % Convert to binary
    encodedDc = [encodedDc, binaryString'];
    % Remove any spaces between the digits (optional)
end

% Loop through each symbol and code
encodedHuff = '';

% Write the Huffman encoding
for idx = 1:length(uniqueValues)
    % Find num and code
    num = dict{idx,1};
    code = dict{idx,2};

    % Calculate the number of bits required to represent the number (signed, sign magnitude)
    if num == 0
        numBits = 2; % If the number is zero, it is represented as having 2 bits
    else
        numBits = floor(log2(abs(num))) + 1 + 1;  % Extra bit for sign magnitude (sign bit)
    end
    % Store the number of bits required as a 4-bit binary number
    numBitsBinary = dec2bin(numBits, 4);  % Store as 4-bit binary
    
    % The number itself in sign magnitude binary representation (with sign bit at the bottom)
    if num < 0
        % Convert negative number to sign magnitude representation
        magnitudeBinary = dec2bin(abs(num), numBits - 1);  % Get the magnitude in binary (without sign bit)
        %disp(magnitudeBinary)
        numberBinary = [magnitudeBinary, '1'];  % Append the sign bit (1) at the bottom (LSB)
    else
        numberBinary = [dec2bin(num, numBits - 1), '0'];  % Positive number, append '0' as sign bit at the bottom
    end
    
    % Number of bits in the Huffman code
    huffmanBits = length(code);
    huffmanBitsBinary = dec2bin(huffmanBits, 6);  % Convert to binary (6-bit format)
    
    % e) The Huffman code itself
    huffmanCode = num2str(code);
    huffmanCode = strrep(huffmanCode, ' ', ''); % Remove spaces

    encodedHuff = [encodedHuff, numBitsBinary, numberBinary, huffmanBitsBinary, huffmanCode];
end

% Binary strings to write
binaryStrings = [encodedHuff, encodedDc, encodedAc];
% Binary strings as rows

% Open a file for writing
filename = sprintf('output/example_a%d_b%d.txt', a, b);

% Open the file with the dynamically generated name
fileID = fopen(filename, 'w'); % 'w' for write mode (overwrites if the file exists)

% Data to write
fprintf(fileID, "%s\n", dec2bin(height, 16));
fprintf(fileID, "%s\n", dec2bin(width, 16));

% Write each row of the array to the file
fprintf(fileID, '%s\n', encodedHuff);
fprintf(fileID, '%s\n', encodedDc);
fprintf(fileID, '%s\n', encodedAc);

% Close the file
fclose(fileID);
    end
end
%         decoded_image = decode_jpeg(encoded_image);
% 
%         error = double(original_image) - double(decoded_image);
%         mse = mean(error(:).^2);
%         rmse = sqrt(mse) / max(original_image(:));
%         RMSE(i, j) = rmse;
% 
%         size_bits = numel(encoded_image) * 8;
%         bpp = size_bits / num_pixels;
%         BPP(i, j) = bpp;
%     end
% 
%     figure;
%     plot(BPP(i, :), RMSE(i, :), '-o', 'LineWidth', 2);
%     xlabel('Bits Per Pixel (BPP)');
%     ylabel('Relative Root Mean Squared Error (RMSE)');
%     title(['RMSE vs BPP for Image ', num2str(i)]);
%     grid on;
% end