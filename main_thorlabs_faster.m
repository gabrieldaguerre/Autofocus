

% ECVFL control using matlab the Optotune Lensdriver 

clear all; 
close all; 
clc; % clear all

% Link tuneable lens
instrreset; %disconnects and deletes all instrument objects
s = serial('COM4'); % depends could change according to the COM fixed
set(s,'BaudRate',115200); % specfication according to Optotune datasheet
fopen(s); % open COM port
fprintf(s,'Start'); % initial command just to check the ECVFL
out = fscanf(s) % prompt print of 'Ready' string in command window

N1=75; %samples 


Contrast_CPP = zeros(1,N1+1);%contrast measurement using the CPP method
Contrast_Teng = zeros(1,N1+1);%contrast measurement using the Teng method
Contrast_Teng_Var = zeros(1,N1+1);%contrast measurement using the Teng-Var method
Contrast_Brenner = zeros(1,N1+1);%contrast measurement using the Brenner method
Contrast_Energy_Lap = zeros(1,N1+1);%contrast measurement using the Energy Lap method
Contrast_Variance_Lap = zeros(1,N1+1);%contrast measurement using the Variance Lap method



% Link camera
% Load TLCamera DotNet assembly. The assembly .dll is assumed to be in the 
% same folder as the scripts.
NET.addAssembly([pwd, '\Thorlabs.TSI.TLCamera.dll']);
disp('Dot NET assembly loaded.');

tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;

% Get serial numbers of connected TLCameras.
serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
disp([num2str(serialNumbers.Count), ' camera was discovered.']);

if (serialNumbers.Count > 0)
    % Open the first TLCamera using the serial number.
    disp('Opening the first camera')
    tlCamera = tlCameraSDK.OpenCamera(serialNumbers.Item(0), false);
    
    % Set exposure time and gain of the camera.
    tlCamera.ExposureTime_us = 20000;
    
    % Check if the camera supports setting "Gain"
    gainRange = tlCamera.GainRange;
    if (gainRange.Maximum > 0)
        tlCamera.Gain = 0;
    end
    
    % Set the FIFO frame buffer size. Default size is 1.
    tlCamera.MaximumNumberOfFramesToQueue = 5;
    
    %figure(1)
    
    % Start continuous image acquisition
    disp('Starting continuous image acquisition.');
    tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
    tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
    tlCamera.Arm;
    tlCamera.IssueSoftwareTrigger;
    maxPixelIntensity = double(2^tlCamera.BitDepth - 1);

    numberOfFramesToAcquire = 1;
    frameCount = 0;
    while frameCount < numberOfFramesToAcquire
        % Check if image buffer has been filled
        if (tlCamera.NumberOfQueuedFrames > 0)
            
            % If data processing in Matlab falls behind camera image
            % acquisition, the FIFO image frame buffer could overflow,
            % which would result in missed frames.
            if (tlCamera.NumberOfQueuedFrames > 1)
                disp(['Data processing falling behind acquisition. ' num2str(tlCamera.NumberOfQueuedFrames) ' remains']);
            end
            
            % For loop : changing the focal power of the tuneable lens
            % (N1 samples)
            for k=1:(N1+1)
                i0=(k-1)*296/N1;%intensiy supplied to thetuneable lens
                current_value =(i0  * 4096) / 293; % from the ECVFL datasheet
                low = 0; 
                high = 0; % variable definition
            % send current value to the optotune lens
            current_value_bin = dec2bin(current_value); % convertion to binary
            size_number = length(current_value_bin);
            while size_number <16% font="">
            current_value_bin = strcat('0',current_value_bin);
            size_number = length(current_value_bin);
            end

            high_bin = current_value_bin(1:8);
            low_bin = current_value_bin(9:end);
            high = bin2dec(high_bin);
            low = bin2dec(low_bin);

            data_command = [hex2dec('41'), hex2dec('77'), high, low, 0, 0];

            crc_value = crc16ibm(data_command, 4); %4=data length
            %   crc16ibm returns a 16 bits number. As before need to split it in 2 8 bit
            % number
            crc_value_bin = dec2bin(crc_value); % convertion to binary
            size_number = length(crc_value_bin);
            while size_number <16% font="">
            crc_value_bin = strcat('0',crc_value_bin);
            size_number = length(crc_value_bin);
            end
            crc_value_high_bin = crc_value_bin(1:8); % upper part
            crc_value_low_bin = crc_value_bin(9:end); % lower part

            crc_value_high = bin2dec(crc_value_high_bin);
            crc_value_low = bin2dec(crc_value_low_bin);

            data_command = [hex2dec('41'), hex2dec('77') high low crc_value_low crc_value_high];
            fwrite(s, data_command, 'uint8');
                        % Get the pending image frame.
            imageFrame = tlCamera.GetPendingFrameOrNull;
            if ~isempty(imageFrame)
                frameCount = frameCount + 1;
                % acquire frame
                % Get the image data as 1D uint16 array
                imageData = uint16(imageFrame.ImageData.ImageData_monoOrBGR);

                disp(['Image frame number: ' num2str(imageFrame.FrameNumber)]);

                % TODO: custom image processing code goes here
                imageHeight = imageFrame.ImageData.Height_pixels;
                imageWidth = imageFrame.ImageData.Width_pixels;
                imageData2D = reshape(imageData, [imageWidth, imageHeight]);
                %I2 = imcrop(im,[x1 y1 x2 y2]) top left bottom right
                I2=imageData2D;
                %I2 = imcrop(imageData2D,[385 674 637 1015]);
                %figure
                %imshow(I2)
                % contrast measurement
                Contrast_CPP(1,k)= Contrast_Measurement_CPP(I2);
                Contrast_Teng(1,k) = Contrast_Measurement_Teng(I2);
                Contrast_Teng_Var(1,k) = Contrast_Measurement_Teng_Var(I2);
                Contrast_Brenner (1,k) = Contrast_Measurement_Brenner(I2);
                Contrast_Energy_Lap (1,k)= Contrast_Measurement_Energy_Laplacian(I2);
                Contrast_Variance_Lap (1,k) = Contrast_Measurement_Variance_Laplacian (I2);
                
                %figure(1),imagesc(imageData2D'), colormap(gray), colorbar
            end
            
            % Release the image frame
            delete(imageFrame);
            end
        end
        drawnow;
    end

% take maximum value
[Max_CPP,index_max_CPP] = max(Contrast_CPP);
[Max_Teng,index_max_Teng] = max(Contrast_Teng);
[Max_Teng_Var,index_max_Teng_Var] = max(Contrast_Teng_Var);
[Max_Brenner,index_max_Brenner] = max(Contrast_Brenner);
[Max_Energy_Lap,index_max_Energy_Lap] = max(Contrast_Energy_Lap);
[Max_Variance_Lap,index_max_Variance_Lap] = max(Contrast_Variance_Lap);

% normalize contrast
Contrast_CPP_norm = Contrast_CPP/ Max_CPP;
Contrast_Teng_norm = Contrast_Teng / Max_Teng;
Contrast_Teng_Var_norm = Contrast_Teng_Var / Max_Teng_Var;
Contrast_Brenner_norm = Contrast_Brenner / Max_Brenner;
Contrast_Energy_Lap_norm = Contrast_Energy_Lap / Max_Energy_Lap;
Contrast_Variance_Lap_norm = Contrast_Variance_Lap / Max_Variance_Lap;

% take current for maximum contrast
i0_CPP = (index_max_CPP-1)*296/N1;
i0_Teng = (index_max_Teng-1)*296/N1;
i0_Teng_Var = (index_max_Teng_Var-1)*296/N1;
i0_Brenner = (index_max_Brenner-1)*296/N1;
i0_Energy_Lap = (index_max_Energy_Lap-1)*296/N1;
i0_Var_Lap = (index_max_Variance_Lap-1)*296/N1;

% compare the results
%t=linspace(1,N1+1,N1+1);
%plot(t,Contrast_CPP_norm,'r',t,Contrast_Teng_norm,'g',t,Contrast_Teng_Var_norm,'b',t,Contrast_Brenner_norm,'m',t,Contrast_Energy_Lap_norm,'c',t,Contrast_Variance_Lap_norm,'k')
%legend('ContrastCPPnorm', 'ContrastTengnorm', 'ContrastTengVarnorm', 'ContrastBrennernorm','ContrastEnergyLapnorm','ContrastVarianceLapnorm')
%xlabel('index ')
%ylabel('Contrast')

% Sub Sampling
G = 20; % size of window between focus

N2=100; % sampling
Contrast_Teng_Var_Sub = zeros(1,N2); %subsampling
Contrast_Brenner_Sub = zeros(1,N2); %subsampling

for k=1:N2
    i0=i0_Teng_Var+(k*G/(N2-2))-(G*N2/(2*(N2-2)));
    current_value =(i0  * 4096) / 293; % from the ECVFL datasheet
    low = 0; 
    high = 0; % variable definition
    % send current to the tunaeble lens
    current_value_bin = dec2bin(current_value); % convertion to binary
    size_number = length(current_value_bin);
    while size_number <16% font="">
        current_value_bin = strcat('0',current_value_bin);
        size_number = length(current_value_bin);
    end

    high_bin = current_value_bin(1:8);
    low_bin = current_value_bin(9:end);
    high = bin2dec(high_bin);
    low = bin2dec(low_bin);

    data_command = [hex2dec('41'), hex2dec('77'), high, low, 0, 0];

    crc_value = crc16ibm(data_command, 4); %4=data length
    %   crc16ibm returns a 16 bits number. As before need to split it in 2 8 bit
    % number
    crc_value_bin = dec2bin(crc_value); % convertion to binary
    size_number = length(crc_value_bin);
    while size_number <16% font="">
        crc_value_bin = strcat('0',crc_value_bin);
        size_number = length(crc_value_bin);
    end
    crc_value_high_bin = crc_value_bin(1:8); % upper part
    crc_value_low_bin = crc_value_bin(9:end); % lower part

    crc_value_high = bin2dec(crc_value_high_bin);
    crc_value_low = bin2dec(crc_value_low_bin);

    data_command = [hex2dec('41'), hex2dec('77') high low crc_value_low crc_value_high];
    fwrite(s, data_command, 'uint8');
                            % Get the pending image frame.
            imageFrame = tlCamera.GetPendingFrameOrNull;
            if ~isempty(imageFrame)
                frameCount = frameCount + 1;
                
                % Get the image data as 1D uint16 array
                imageData = uint16(imageFrame.ImageData.ImageData_monoOrBGR);

                disp(['Image frame number: ' num2str(imageFrame.FrameNumber)]);

                % TODO: custom image processing code goes here
                imageHeight = imageFrame.ImageData.Height_pixels;
                imageWidth = imageFrame.ImageData.Width_pixels;
                imageData2D = reshape(imageData, [imageWidth, imageHeight]);
                % calculate contrast in the image center
                I2 = imcrop(imageData2D,[385 674 637 1015]);
                %I2=imageData2D;
                %figure
                %imshow(I2,'DisplayRange',[0  255])
                Contrast_Teng_Var_Sub(1,k) = Contrast_Measurement_Teng_Var(I2);
                Contrast_Brenner_Sub(1,k) = Contrast_Measurement_Brenner(I2);
                %figure(1),imagesc(imageData2D'), colormap(gray), colorbar
            end
            
            % Release the image frame
            delete(imageFrame);
 end
 end
 drawnow;
% take max    
[Max_Teng_Var_Sub,index_max_Teng_Var_Sub] = max(Contrast_Teng_Var_Sub);
[Max_Brenner_Sub,index_max_Brenner_Sub] = max(Contrast_Brenner_Sub);
% take current for maximum contrast
i0_Teng_Var_Sub = i0_Teng_Var+(index_max_Teng_Var_Sub*G/(N2-2))-(G*N2/(2*(N2-2)));
i0_Brenner_Sub = i0_Brenner+(index_max_Brenner_Sub*G/(N2-2))-(G*N2/(2*(N2-2)));


    


% Set to Best focus position
    i0=i0_Teng_Var_Sub;
    current_value =(i0  * 4096) / 293; % from the ECVFL datasheet
    low = 0; 
    high = 0; % variable definition
    current_value_bin = dec2bin(current_value); % convertion to binary
    size_number = length(current_value_bin);
    while size_number <16% font="">
        current_value_bin = strcat('0',current_value_bin);
        size_number = length(current_value_bin);
    end

    high_bin = current_value_bin(1:8);
    low_bin = current_value_bin(9:end);
    high = bin2dec(high_bin);
    low = bin2dec(low_bin);

    data_command = [hex2dec('41'), hex2dec('77'), high, low, 0, 0];

    crc_value = crc16ibm(data_command, 4); %4=data length
    %   crc16ibm returns a 16 bits number. As before need to split it in 2 8 bit
    % number
    crc_value_bin = dec2bin(crc_value); % convertion to binary
    size_number = length(crc_value_bin);
    while size_number <16% font="">
        crc_value_bin = strcat('0',crc_value_bin);
        size_number = length(crc_value_bin);
    end
    crc_value_high_bin = crc_value_bin(1:8); % upper part
    crc_value_low_bin = crc_value_bin(9:end); % lower part

    crc_value_high = bin2dec(crc_value_high_bin);
    crc_value_low = bin2dec(crc_value_low_bin);

    data_command = [hex2dec('41'), hex2dec('77') high low crc_value_low crc_value_high];
    fwrite(s, data_command, 'uint8');
                                % Get the pending image frame.
            imageFrame = tlCamera.GetPendingFrameOrNull;
            if ~isempty(imageFrame)
                frameCount = frameCount + 1;
                
                % Get the image data as 1D uint16 array
                imageData = uint16(imageFrame.ImageData.ImageData_monoOrBGR);

                disp(['Image frame number: ' num2str(imageFrame.FrameNumber)]);

                % TODO: custom image processing code goes here
                imageHeight = imageFrame.ImageData.Height_pixels;
                imageWidth = imageFrame.ImageData.Width_pixels;
                imageData2D = reshape(imageData, [imageWidth, imageHeight]);

                %figure(1),imagesc(imageData2D'), colormap(gray), colorbar
            end
            
            % Release the image frame
            delete(imageFrame);


    % Stop continuous image acquisition
    disp('Stopping continuous image acquisition.');
    tlCamera.Disarm;
    
    % Release the TLCamera
    disp('Releasing the camera');
    tlCamera.Dispose;
    delete(tlCamera);


% Release the serial numbers
delete(serialNumbers);

% Release the TLCameraSDK.
tlCameraSDK.Dispose;
delete(tlCameraSDK);


instrreset;

