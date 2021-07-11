clear;
close all;
clc;
profile on
UpFolder = fileparts(pwd);
addpath(fullfile(UpFolder, 'libraries/l1magic/Optimization'));
addpath(fullfile(UpFolder, 'libraries/l1magich/Measurements'));
addpath(fullfile(UpFolder, 'libraries/l1magic/Data'));
addpath(fullfile(UpFolder, 'libraries/spgl1-1.9'));
addpath(fullfile(UpFolder, 'libraries/SL0'));
addpath(fullfile(UpFolder, 'libraries/TwIST_v2'));
addpath(fullfile(UpFolder, 'libraries/lzw'));
addpath(fullfile(UpFolder, 'libraries/pureAC'));
addpath(fullfile(UpFolder, 'libraries/'));
addpath(fullfile(UpFolder, 'libraries/BallLabsAlgo'));
addpath(fullfile(UpFolder, 'Frame\'));

imageOriginalPath = 'C:\Users\Jay\Desktop\Research\Legacy_IEEE_ACCESS_CS_based_image_compression';
imageFiles = [dir(fullfile(imageOriginalPath,'*png'));
              dir(fullfile(imageOriginalPath,'*tiff'));
              dir(fullfile(imageOriginalPath,'*tif'));
              dir(fullfile(imageOriginalPath,'*jpg'));
              dir(fullfile(imageOriginalPath,'*bmp'));
              dir(fullfile(imageOriginalPath,'*mat'))];
numFiles = length(imageFiles);
%___SIMULATION SETUPS___
simulation_parameter.macro_block_size                = 8;
simulation_parameter.n                               = simulation_parameter.macro_block_size^2; % NOTE: small error still present after increasing m;
simulation_parameter.measurement_matrix_lists        = [48 32 16 8];
simulation_parameter.bit                             = [3 4 5 6];
simulation_parameter.measurement_matrix_construction = 'binary_random';
simulation_parameter.reconstruction_algorithm        = 'l1_eq_pd';
simulation_parameter.transformation_algorithm        = 'idct';
simulation_parameter.color_mode                      = 'a_grey';

for matrix_loop = 1:length(simulation_parameter.measurement_matrix_lists)
    for quantization_bit = 1:length(simulation_parameter.bit)
        reset = 1;
        switch simulation_parameter.measurement_matrix_lists(matrix_loop)
            case {64}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 1;
            case {48}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.75;
            case {32}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.50;
            case {16}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.25;        
            case {8}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.1250;
            case {4}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.0625;
            case {2}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.03125;
            case {1}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.015625;
            case {0.5}
                simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
                simulation_parameter.sampling_rate = 0.0078125;
        end
        
        switch simulation_parameter.measurement_matrix_construction
            case 'binary_random'
                temp_matrix = load('random_matrix.mat');
                simulation_parameter.phi = temp_matrix.ans(1:simulation_parameter.m,1:simulation_parameter.n);
                %simulation_parameter.phi = randi([0, 1], [simulation_parameter.m, simulation_parameter.n]); % This will give you a friendly measurement matrix (M must smaller than N)
            case 'binary_walsh_hadamard'
                hadamard_matrix          = hadamard(simulation_parameter.n);
                HadIdx                   = 0:simulation_parameter.n-1;       % Hadamard index
                M                        = log2(simulation_parameter.n)+1;    % Number of bits to represent the index
                binHadIdx                = fliplr(dec2bin(HadIdx,M))-'0';     % Bit reversing of the binary index
                binSeqIdx                = zeros(simulation_parameter.n,M-1); % Pre-allocate memory
                for k = M:-1:2
                    % Binary sequency index
                    binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
                end
                SeqIdx                   = binSeqIdx*pow2((M-1:-1:0)');   % Binary to integer sequency index
                walshMatrix              = hadamard_matrix(SeqIdx+1,:);   % 1-based indexing
                simulation_parameter.phi = max(walshMatrix(1:simulation_parameter.m,1:simulation_parameter.n), 0);
            case 'binary_walsh_hadamard_zigzag'
                hadamard_matrix          = hadamard(simulation_parameter.n);
                HadIdx                   = 0:simulation_parameter.n-1;        % Hadamard index
                M                        = log2(simulation_parameter.n)+1;    % Number of bits to represent the index
                binHadIdx                = fliplr(dec2bin(HadIdx,M))-'0';     % Bit reversing of the binary index
                binSeqIdx                = zeros(simulation_parameter.n,M-1); % Pre-allocate memory
                for k = M:-1:2
                    % Binary sequency index
                    binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
                end
                SeqIdx                   = binSeqIdx*pow2((M-1:-1:0)');   % Binary to integer sequency index
                walshMatrix              = hadamard_matrix(SeqIdx+1,:);   % 1-based indexing
                simulation_parameter.phi = max(walshMatrix,0);
                sub_phi.x                = zeros(1, size(simulation_parameter.phi,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
                sub_phi.y                = zeros(1, size(simulation_parameter.phi,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
                simulation_parameter.phi = mat2cell(simulation_parameter.phi, sub_phi.x, sub_phi.y);
                %___Re-order to AC DC
                t=0;
                l=size(simulation_parameter.phi);
                sum_=l(2)*l(1);  %calculating the M*N
                for d=2:sum_
                 c=rem(d,2);  %checking whether even or odd
                    for i=1:l(1)
                        for j=1:l(2)
                            if((i+j)==d)
                                t=t+1;
                                if(c==0)
                                simulation_parameter.phi_temp(t,:) = simulation_parameter.phi{j,d-j}(:)';
                                else          
                                simulation_parameter.phi_temp(t,:) = simulation_parameter.phi{d-j,j}(:)';
                                end
                             end    
                         end
                     end
                end
                simulation_parameter.phi = max(simulation_parameter.phi_temp(1:simulation_parameter.m,1:simulation_parameter.n),0);
            case 'binary_walsh_hadamard_hubert'
            case 'binary_hadamard'
                hadamard_matrix          = hadamard(simulation_parameter.n);  
                simulation_parameter.phi = max(hadamard_matrix,0);
                sub_phi.x                = zeros(1, size(simulation_parameter.phi,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
                sub_phi.y                = zeros(1, size(simulation_parameter.phi,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
                simulation_parameter.phi = mat2cell(simulation_parameter.phi, sub_phi.x, sub_phi.y);
                %___Re-order to AC DC
                t=0;
                l=size(simulation_parameter.phi);
                sum_=l(2)*l(1);  %calculating the M*N
                for d=2:sum_
                 c=rem(d,2);  %checking whether even or odd
                    for i=1:l(1)
                        for j=1:l(2)
                            if((i+j)==d)
                                t=t+1;
                                if(c==0)
                                    simulation_parameter.phi_temp(t,:) = simulation_parameter.phi{j,d-j}(:)';
                                else          
                                    simulation_parameter.phi_temp(t,:) = simulation_parameter.phi{d-j,j}(:)';
                                end
                             end    
                         end
                     end
                end
                simulation_parameter.phi = max(simulation_parameter.phi_temp(1:simulation_parameter.m,1:simulation_parameter.n),0);
        end

        mat = simulation_parameter.phi;  % Your sample matrix
        [r, c] = size(mat);                          % Get the matrix size
        imagesc((1:c)+0.5, (1:r)+0.5, mat);          % Plot the image
        colormap(gray);                              % Use a gray colormap
        axis equal                                   % Make axes grid sizes equal
        set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
                 'XLim', [1 c+1], 'YLim', [1 r+1], ...
                 'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
        %___THETA___
        %___NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues___
        simulation_parameter.theta = zeros(simulation_parameter.m,simulation_parameter.n);
        for theta_loop = 1:simulation_parameter.n
            ek = zeros(1,simulation_parameter.n);
            ek(theta_loop) = 1;
            switch simulation_parameter.transformation_algorithm
                case 'idct'
                    simulation_parameter.psi = idct2(ek)';
                case 'ifwht'
                    simulation_parameter.psi = ifwht(ek)';
            end
            simulation_parameter.theta(:,theta_loop) = simulation_parameter.phi*simulation_parameter.psi;
        end
            frame_index = 1;
            for frame_number = 1:1
                frame_number
                %___LOAD IMAGE___
                load_frame = imread(imageFiles(frame_number).name);
                if(strcmp(simulation_parameter.color_mode,'rgb') || strcmp(simulation_parameter.color_mode,'RGB'))
                    frame = load_frame(:,:,:);
                    frame = frame(1:640, 1:640, :);
                    plane = 3;
                elseif(strcmp(simulation_parameter.color_mode,'gray') || strcmp(simulation_parameter.color_mode,'GRAY'))
                    frame = rgb2gray(load_frame);
                    frame = frame(1:600, 1:800, :);
                    plane = 1;
                else
                    frame = load_frame;
                    %frame = frame(1:640, 1:640, :);
                    plane = 1;
                end

                %___RESET STATE___
                sub_block.x = zeros(1, size(frame,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
                sub_block.y = zeros(1, size(frame,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
                for i = 1:plane
                    frame_temp(:,:,i) = mat2cell(frame(:,:,i), sub_block.x, sub_block.y);
                end
                
                %___RESET STATE___
                if(reset == 1)
                    reset                               = 0;
                    bits                                = 0;
                    quantization.quantization_parameter = zeros(simulation_parameter.m, 1);
                    y.up_encoder                        = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                    y.up_encoder_SDPC                   = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                    y.up_encoder_JB                     = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                    y.up_encoder_new                    = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                    y.up_decoder                        = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                    y.left_encoder                      = zeros(simulation_parameter.m, 1);
                    y.left_encoder_SDPC                 = zeros(simulation_parameter.m, 1);
                    y.left_encoder_JB                   = zeros(simulation_parameter.m, 1);
                    y.left_encoder_new                  = zeros(simulation_parameter.m, 1);
                    y.up_left_encoder                   = zeros(simulation_parameter.m, 1);
                    y.up_left_encoder_SDPC              = zeros(simulation_parameter.m, 1);
                    y.up_left_encoder_JB                = zeros(simulation_parameter.m, 1);
                    y.up_left_encoder_new               = zeros(simulation_parameter.m, 1);
                    y.up_right_encoder                  = zeros(simulation_parameter.m, 1);
                    y.dc_encoder                        = zeros(simulation_parameter.m, 1);
                    y.left_decoder                      = zeros(simulation_parameter.m, 1);
                    y.up_left_decoder                   = zeros(simulation_parameter.m, 1);
                    y.up_right_decoder                  = zeros(simulation_parameter.m, 1);
                    y.dc_decoder                        = zeros(simulation_parameter.m, 1);
                    y.predicted_encoder                 = zeros(simulation_parameter.m, 1);
                    y.measurement                       = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    y.residual                          = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    y.combine                           = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    y.quantized                         = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    y.dequantized                       = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    temp_padding                        = zeros(size(frame,1)+2, size(frame,2)+2, plane);
                    res_temp_padding                    = zeros(size(frame,1)+2, size(frame,2)+2, plane);
                    reconstructed_image                 = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    res_reconstructed_image             = cell(size(sub_block.x,2), size(sub_block.y,2), plane);
                    modes                               = zeros(size(frame,1)/simulation_parameter.macro_block_size, size(frame,2)/simulation_parameter.macro_block_size, plane);
                    res_video_buffer                    = zeros(size(frame,1), size(frame,2), plane, 100);
                    video_buffer                        = zeros(size(frame,1), size(frame,2), plane, 100);
                    for i = 1:size(sub_block.x, 2)
                        for j = 1:size(sub_block.y, 2)
                            for k = 1:plane
                                y.measurement{i,j,k}           = zeros(simulation_parameter.m, 1);
                                y.combine{i,j,k}               = zeros(simulation_parameter.m, 1);
                                y.residual{i,j,k}              = zeros(simulation_parameter.m, 1);
                                y.quantized{i,j,k}             = zeros(simulation_parameter.m, 1);
                                y.dequantized{i,j,k}           = zeros(simulation_parameter.m, 1);
                                reconstructed_image{i,j,k}     = zeros(simulation_parameter.m, 1);
                                res_reconstructed_image{i,j,k} = zeros(simulation_parameter.m, 1);
                            end
                        end
                    end
                end

                %___THE RANDOM PROJECTION___
                disp('Random prejection...');
                for k = 1:plane
                    for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                        for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                           one_block_image(:,:,k)         = reshape(frame_temp{i,j,k},simulation_parameter.macro_block_size^2,1);
                           y.measurement{i,j,k}           = BCS_encoder(double(one_block_image(:,:,k)), simulation_parameter.phi); %___Sampling
                        end
                    end
                end

                residual_entropy_SQ = 0;
                residual_entropy_DPCM = 0;
                residual_entropy_SDPC = 0;
                residual_entropy_JB = 0;
                residual_entropy_new = 0;

                disp('Intra-Frame Coding...');
                residual_entropy = 0;
                previous = zeros(simulation_parameter.m, 1);
                for k = 1:plane
                    for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                        for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                            y.residual{i,j,k}                               = y.measurement{i,j,k}; 
                            %___SQ coding___
                            bit_depth                                       = simulation_parameter.bit(quantization_bit); % compute number of levles    
                            quantization.quantization_parameter_SQ(i,j,k)   = (max(y.residual{i,j,k})-min(y.residual{i,j,k}))/2^bit_depth; % compute step size of uniform quantizer
                            y.quantized_SQ{i,j,k}                           = round(y.residual{i,j,k}/quantization.quantization_parameter_SQ(i,j,k));  %compute quantized image
                            y.quantized_SQ{i,j,k}(isinf(y.quantized_SQ{i,j,k})|isnan(y.quantized_SQ{i,j,k})) = 0; % Replace NaNs and infinite values with zeros
                            residual_entropy_SQ                             = residual_entropy_SQ + Measurement_Entropy(y.quantized_SQ{i,j,k}, size(frame,1)*size(frame,2));
                            y.dequantized_SQ{i,j,k}                         = (y.quantized_SQ{i,j,k}*quantization.quantization_parameter_SQ(i,j,k));
                            previous                                        = y.dequantized_SQ{i,j,k};
                            
                            x_hat_SQ{i,j,k}                                 = BCS_reconstruction(y.dequantized_SQ{i,j,k}, ...
                                                                                                 simulation_parameter.theta, ...
                                                                                                 simulation_parameter.phi, ...
                                                                                                 simulation_parameter.transformation_algorithm , ...
                                                                                                 simulation_parameter.reconstruction_algorithm, ...
                                                                                                 simulation_parameter.macro_block_size);  
                        end
                    end
                end
                
                disp('Finished...');
                %___QUATITATIVE MATRICES___
                quantitative_metrices.psnr_SQ(matrix_loop, quantization_bit)             = PSNR(uint8(cell2mat(x_hat_SQ(:,:,:,frame_number))), frame);
                quantitative_metrices.ssim_SQ(matrix_loop, quantization_bit)             = ssim(uint8(cell2mat(x_hat_SQ(:,:,:,frame_number))), frame);
                quantitative_metrices.entropy_bits_SQ(matrix_loop, quantization_bit)     = residual_entropy_SQ;
            end
    end
end
colorstring = 'kbgry';
figure;
for i = 1:length(simulation_parameter.bit) 
    plot(quantitative_metrices.entropy_bits_SQ(:,i), quantitative_metrices.ssim_SQ(:,i), 'Color', colorstring(i), 'LineWidth', 2);
    hold on
end
title('Parrots');
xlabel('Bitrate (bpp)', 'Fontsize', 10);
ylabel('SSIM', 'Fontsize', 10);
legend('75%', '50%', '25%', '12.5%','Fontsize', 10,'Location', 'SouthEast');

profile report
profile off