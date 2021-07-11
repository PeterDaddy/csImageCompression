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
addpath(fullfile(UpFolder, 'libraries'));
addpath(fullfile(UpFolder, 'Frame/'));

imageOriginalPath = 'E:\Research\Frame/';
imageFiles = [dir(fullfile(imageOriginalPath,'*png'));
              dir(fullfile(imageOriginalPath,'*tiff'));
              dir(fullfile(imageOriginalPath,'*tif'));
              dir(fullfile(imageOriginalPath,'*jpg'));
              dir(fullfile(imageOriginalPath,'*bmp'));
              dir(fullfile(imageOriginalPath,'*mat'))];
numFiles = length(imageFiles);
%___SIMULATION SETUPS___
simulation_parameter.macro_block_size                = 16;
simulation_parameter.n                               = simulation_parameter.macro_block_size^2; % NOTE: small error still present after increasing m;
simulation_parameter.measurement_matrix_lists        = [128];
simulation_parameter.measurement_matrix_construction = 'binary_walsh_hadamard';
simulation_parameter.reconstruction_algorithm        = 'l1_eq_pd';
simulation_parameter.transformation_algorithm        = 'ifwht';
simulation_parameter.color_mode                      = 'gray';

for matrix_loop = 1:length(simulation_parameter.measurement_matrix_lists)
    reset = 1;
    switch simulation_parameter.measurement_matrix_lists(matrix_loop)
        case {256}
            simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
            simulation_parameter.sampling_rate = 1;
        case {192}
            simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
            simulation_parameter.sampling_rate = 0.75;
        case {128}
            simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
            simulation_parameter.sampling_rate = 0.50;
        case {64}
            simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
            simulation_parameter.sampling_rate = 0.25;        
        case {32}
            simulation_parameter.m             = simulation_parameter.measurement_matrix_lists(matrix_loop);
            simulation_parameter.sampling_rate = 0.1250;
    end

    switch simulation_parameter.measurement_matrix_construction
        case 'binary_random'
            simulation_parameter.phi = randi([0, 1], [simulation_parameter.m, simulation_parameter.n]); % This will give you a friendly measurement matrix (M must smaller than N)
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
        case 'binary_hadamard'
            hadamard_matrix          = hadamard(simulation_parameter.n);
            simulation_parameter.phi = max(hadamard_matrix(1:simulation_parameter.m,1:simulation_parameter.n), 0);
    end

    %___THETA___
    %___NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues___
    theta = zeros(simulation_parameter.m,simulation_parameter.n);
    for theta_loop = 1:simulation_parameter.n
        ek = zeros(1,simulation_parameter.n);
        ek(theta_loop) = 1;
        switch simulation_parameter.transformation_algorithm
            case 'idct'
                simulation_parameter.psi = idct2(ek)';
            case 'ifwht'
                simulation_parameter.psi = ifwht(ek, length(ek),'sequency')';
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
                frame = imresize(frame,[720, 1280]);
                plane = 3;
            else
                frame = rgb2gray(load_frame);
                %frame = imresize(frame,[720, 1280]);
                plane = 1;
            end

            %___RESET STATE___
            sub_block.x = zeros(1, size(frame,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
            sub_block.y = zeros(1, size(frame,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;

            for i = 1:plane
                frame_temp(:,:,i) = mat2cell(im2double(frame(:,:,i)), sub_block.x, sub_block.y);
            end
            %___RESET STATE___
            if(reset == 1)
                reset                               = 0;
                bits                                = 0;
                quantization.quantization_parameter = zeros(simulation_parameter.m, 1);
                y.up_encoder                        = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                y.up_decoder                        = zeros((simulation_parameter.m), size(frame,2)/simulation_parameter.macro_block_size);
                y.left_encoder                      = zeros(simulation_parameter.m, 1);
                y.up_left_encoder                   = zeros(simulation_parameter.m, 1);
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
                quantitative_metrices.coded_bit     = zeros(1,3);
                quantitative_metrices.psnr          = zeros(1,3);
                quantitative_metrices.ssim          = zeros(1,3);
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
    
            dec = load('cs_compressed.cseg2020', '-mat'); % Read compress data from the file
            ac_r.message_coded_1  = dec.Header(1).message_coded_1;
            ac_r.message_counts_1 = dec.Header(2).message_counts_1;
            ac_r.message_table_1  = dec.Header(3).message_table_1;
            ac_r.message_coded_2  = dec.Header(4).message_coded_2;
            ac_r.message_counts_2 = dec.Header(5).message_counts_2;
            ac_r.message_table_2  = dec.Header(6).message_table_2;

            ac_r.modes_coded_1    = dec.Header(7).modes_coded_1;
            ac_r.modes_counts_1   = dec.Header(8).modes_counts_1;
            ac_r.modes_table_1    = dec.Header(9).modes_table_1;
            ac_r.modes_coded_2    = dec.Header(10).modes_coded_2;
            ac_r.modes_counts_2   = dec.Header(11).modes_counts_2;
            ac_r.modes_table_2    = dec.Header(12).modes_table_2;
            
            %----------------Reciever part--------------
            disp('Arithmatic Decoding Started');
            ac_r.rle_message{1} = Arith_Decode(ac_r.message_coded_1,ac_r.message_counts_1,ac_r.message_table_1);
            ac_r.rle_message{2} = Arith_Decode(ac_r.message_coded_2,ac_r.message_counts_2,ac_r.message_table_2);
            
            ac_r.rle_modes{1}   = Arith_Decode(ac_r.modes_coded_1,ac_r.modes_counts_1,ac_r.modes_table_1);
            ac_r.rle_modes{2}   = Arith_Decode(ac_r.modes_coded_2,ac_r.modes_counts_2,ac_r.modes_table_2);
            
            ac_r.dec_message    = rle(ac_r.rle_message);
            ac_r.dec_mode       = rle(ac_r.rle_modes);
            ac_r.modes          = reshape(ac_r.dec_mode, [32,32]);
            disp(['Recieved message is    ',ac_r.dec_message]);

            %___Package translation___
            disp('Package translation to decoding format');
            first_member_position = 1;
            last_member_position = simulation_parameter.m;
            for k = 1:plane
                for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                    for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                        y.streaming_received{i,j,k} = ac_r.dec_message(first_member_position:last_member_position)';
                        first_member_position       = last_member_position + 1;
                        last_member_position        = last_member_position + simulation_parameter.m;        
                    end
                end
            end
                        
            for k = 1:plane
                for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                    for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                       %___DE-QUANTIZATION___
                       y.dequantized{i,j,k} = round(y.streaming_received{i,j,k}.*256);
                       if(ac_r.modes(i,j,k) == 1)
                          y.combine{i,j,k} = y.dequantized{i,j,k} + y.left_decoder;
                       elseif(ac_r.modes(i,j,k) == 2)
                          y.combine{i,j,k} = y.dequantized{i,j,k} + y.up_decoder(:,j);
                       elseif(ac_r.modes(i,j,k) == 3)
                          y.combine{i,j,k} = y.dequantized{i,j,k} + y.up_left_decoder;
                       elseif(ac_r.modes(i,j,k) == 4)
                          y.combine{i,j,k} = y.dequantized{i,j,k} + y.up_right_decoder;
                       elseif(ac_r.modes(i,j,k) == 5)
                          y.combine{i,j,k} = y.dequantized{i,j,k} + round(y.left_decoder + y.up_decoder(:,j))/2;   
                       else
                          y.combine{i,j,k} = y.dequantized{i,j,k};
                       end

                       %__Next iteration__
                       y.left_decoder         = y.combine{i,j,k};
                       y.up_decoder(:,j)      = y.combine{i,j,k};
                       if((i>1)&&(j>2))
                           y.up_left_decoder  = y.up_decoder(:,j-1);
                       else
                           y.up_left_decoder  = zeros(simulation_parameter.m,1);
                       end
                       if((i>1)&&(j<size(frame,2)/simulation_parameter.macro_block_size))
                           y.up_right_decoder = y.up_decoder(:,j+1);
                       else
                           y.up_right_encoder = zeros(simulation_parameter.m,1);
                       end

                      reconstructed_image{i,j,k} = BCS_reconstruction(y.measurement{i,j,k}, ...
                                                                      simulation_parameter.theta, ...
                                                                      simulation_parameter.phi, ...
                                                                      simulation_parameter.transformation_algorithm, ...
                                                                      simulation_parameter.reconstruction_algorithm, ...
                                                                      simulation_parameter.macro_block_size);
                    end
                end
            end
%             %___FINAL PROCESS ZONE___%
%             for k = 1:plane
%                 if(plane == 3)
%                    temp_padding(:,:,k) = padarray(floor(cell2mat(reconstructed_image(:,:,k))),[1 1],'symmetric','both');
%                 else
%                    temp_padding = padarray(floor(cell2mat(reconstructed_image(:,:))),[1 1],'symmetric','both'); %padding
%                 end
% 
%                 %___Overlapped Filtering___ Y X
%                 for i = 2:size(temp_padding,1)-1
%                     for j = 2:size(temp_padding,2)-1
%                         video_buffer(i-1,j-1,k,frame_number) = floor((temp_padding(i,j,k)+temp_padding(i,j-1,k)+temp_padding(i,j+1,k))/3);
%                     end
%                 end
%             end
        %___WINER___
        video_buffer(:,:,:,frame_number) = wiener2(cell2mat(reconstructed_image(:,:)),[3 3]);

        %___RESET FOR NEW FRAME___
%         y.buffer_up_encoder      = zeros((m), size(frame,2)/macro_block_size);
%         y.buffer_left_encoder    = zeros(m, 1);
%         y.buffer_dc_encoder      = zeros(m, 1);
%         y.buffer_cp_encoder      = (zeros(m, 1));
%         y.buffer_up_decoder      = zeros((m), size(frame,2)/macro_block_size);
%         y.buffer_left_decoder    = zeros(m, 1);
%         y.buffer_dc_decoder      = zeros(m, 1);
%         y.buffer_cp_decoder      = (zeros(m, 1));
%         y.buffer_cp_encoder(1)   = 0;
%         y.buffer_cp_encoder(2:m) = 0;
%         y.buffer_cp_decoder(1)   = 0;
%         y.buffer_cp_decoder(2:m) = 0;
%         modes                    = zeros(size(frame,1)/macro_block_size, size(frame,2)/macro_block_size, plane);

        %___QUATITATIVE MATRICES___
%         quantitative_metrices.storage_size(frame_number) = ByteSize(Header);
         quantitative_metrices.psnr(frame_number)         = PSNR(video_buffer(:,:,:,frame_number)*255, im2double(frame)*255);
         quantitative_metrices.ssim(frame_number)         = ssim(video_buffer(:,:,:,frame_number), im2double(frame));
         bits                                             = 0;
      end

%     save(fullfile(strcat('PSNR_RedKayak_', num2str(sampling_rate),'_Qp_', num2str(2^log2(M)), '.mat')), 'quantitative_metrices.psnr');
%     save(fullfile(strcat('SSIM_RedKayak_', num2str(sampling_rate),'_Qp_', num2str(2^log2(M)), '.mat')), 'quantitative_metrices.ssim');
%     save(fullfile(strcat('BPP_RedKayak_' , num2str(sampling_rate),'_Qp_', num2str(2^log2(M)), '.mat')), 'quantitative_metrices.bpp');
% 
%     video_out = VideoWriter(fullfile(strcat('RedKayak', ...
%                                             '_Frame_Skip_', num2str(sampling_rate), ...
%                                             '_', measurement_matrix_construction, ...
%                                             '_', quantitative_metrices.reconstruction_algorithm, ...
%                                             '_', quantitative_metrices.transformation_algorithm, ...
%                                             '_', color_mode, ...
%                                             '_Linear_Filter', ...
%                                             '.avi')), 'Uncompressed AVI'); %create the video object
%     video_out.FrameRate = 30;                                    
%     open(video_out); %open the file for writing
%     for loop = 1:frame_number-1
%        writeVideo(video_out, uint8(video_buffer(:,:,:,loop)));
%     end
%     close(video_out);
end
profile report
profile off