%------------------------------------------------------
%   2021 RT516 - Medical Imaging
%   Term Project 2, Part 2 - Ultrasound Imaging
%   Part II : Reconstruct raw RF data into an image
%   202123008 KIM Jinmin (M.S. Candidate)
%   Department of Robotics Engineering
%------------------------------------------------------
clear all
close all
clc

%% 이미지 로드, 변수 지정

%Read in the image data
[imagedata, numVectors, numElements, numSamples] = readBinData('imageData.bin');

% Determine sizes of image data (depth samples x number of elements (channels) x number of transmit lines)
[samplenum, elemnum, Alinenum] = size(imagedata);

receivenum=10;  % number of receive foci in depth for dynamic receive focusing- this can be changed but does not need to be
times = zeros(1,samplenum-9);   % make an array for storing delay times

samplestart = zeros(1,receivenum);  % make an array for storing start sample number for a given receive focal zone in depth
sampleend = zeros(1,receivenum);    % make an array for storing end sample number for a given receive focal zone in depth
z = zeros(1,receivenum);            % depth of each receive focal zone in mm
delay = zeros(Alinenum, elemnum, receivenum); % matrix to store computed delay values
beamspacing=0.177;                  % separation of Rx beams in mm
txspacing=0.201;                    % separation of Tx beams in mm
center = 96*.201;                   % distance across half of the array
c = 1540000;                        % speed of sound in mm/s - You will change this in Question 2
total_depth = 37;                   % max depth to reconstruct in mm
fs = 40e6;                          % Sampling frequency (Hz) (Axial axis)

Xf=-beamspacing*20:beamspacing:beamspacing*20;  % lateral (x) locations to reconstruct - 41 total
Xi=-txspacing*(elemnum/2-.5):txspacing:txspacing*(elemnum/2-.5);    % lateral (x) locations of transducer elements  - 192 total
Zi=zeros(1,192);    % depth (z) locations of transducer elements

%%
for k = 1:receivenum %for the number of receive focal zones specified
    % ----------------------------------------------------
    % Focal zone calculation
    % ----------------------------------------------------    
    depth_increment = total_depth/((receivenum+1));%the distance between focal zones in depth
    z(k) = depth_increment*(k);
    if k == 1
        samplestart(k) = 1;%for the first focal zone, start with sample 1
    else
        samplestart(k) = round(fs*2*(z(k) - depth_increment/2)/c);% for other focal zones, start at a depth determined by the focal zone size and sampling frequency, which is 40 MHz
    end
    sampleend(k) = round(fs*2*(z(k)+depth_increment/2)/c);%compute end of each focal zone

    % Receive focal point calculation
    Rf=sqrt(Xf.^2 + z(k).^2); % distance between focal point and phase center (i.e. (0,0) on transducer)
    for a_line = 1:41
        for element = 1:elemnum
            %---------------------------------------------------------------------------------
            %YOU NEED TO ADD CODE HERE TO COMPUTE THE DELAY VALUES!
            
            temp = element * a_line * k / c / fs;  % temp is the delay value. I have to write the code to calculate focal delay.
     
            %---------------------------------------------------------------------------------
            delay(a_line,element,k) = temp; %compute focal delay for this line, element, and depth
        end
        %---------------------------------------------------------------------------------
        %IT IS SUGGESTED TO ADD SOME CODE HERE TO SUBTRACT THE MINIMUM
        %DELAY VALUE ACROSS ALL ELEMENTS, BUT ULTIMATELY YOU CAN DO THIS HOWEVER
        %YOU WOULD LIKE
        %---------------------------------------------------------------------------------
        % delay(a_line,:,k) = 0;% subtract minimum delay.  This structure stores the delay value in sec
    end
end
delaynum = round(delay./(1/fs));%This structure stores the delay value in samples
zsample = round((fs*2*z/c));%This stores the sample number associated with each Rx focal zone

%%
%-------------------------------------------------------
% IN PART 3: Apply your apodization window HERE - it will need to be applied to ALL
% the data in imagedata, i.e. your window function should be a vector of
% length 192 (for 192 elements) and should be applied to all 41 transmit
% lines and all depth samples.
% For Parts 1 and 2, you do need to add any code here and the image will
% still be reconstructed successfully.
%-------------------------------------------------------

%-------------------------------------------------------
final_image = zeros(samplenum-9-max(max(delaynum(:,:,receivenum))),Alinenum,receivenum);%Create storage for final image (i.e. after beamforming)
for Aline = 1:Alinenum
    for k = 1:receivenum
        for sample = 10:(samplenum - max(max(delaynum(:,:,k))))
            time = (sample - 10)*(1/fs);
            times(sample-9) = time;
            depth_mm = time*c/2;
            for element = 1:elemnum
                %---------------------------------------------------------------------------------
                %YOU NEED TO ADD CODE HERE TO DO SUMMATION
                
                Summ = 0; %this is where beamforming (summation part) is applied
                
                %---------------------------------------------------------------------------------
                final_image(sample-9,Aline,k) = Summ; %this is where beamforming (summation part) is applied
            end
        end
    end
end
sampleend(receivenum) = samplenum-9-max(max(delaynum(:,:,receivenum)));%compute
final_dynimage = zeros(samplenum-9-max(max(delaynum(:,:,receivenum))),Alinenum);
for k = 1:receivenum
    if k == 1
        final_dynimage(samplestart(k):sampleend(k),:) = final_image(samplestart(k):sampleend(k),:,k);
    else
        final_dynimage(samplestart(k)-1:sampleend(k),:) = final_image(samplestart(k)-1:sampleend(k),:,k);
    end
end
final_image_hilbert = abs(hilbert(final_dynimage));
final_image_normalized = 20*log10(final_image_hilbert./max(final_image_hilbert(:)));
final_image_compression = 20*log10(final_image_hilbert./max(final_image_hilbert(:)));
final_image_compression(zsample(1)-2:zsample(1)+2,:) = 10000;
final_image_compression(samplestart(1):samplestart(1)+15,:) = 500;
for k = 2:receivenum
    final_image_compression(zsample(k)-2:zsample(k)+2,:) = 10000;
    if k == receivenum
        final_image_compression(sampleend(k)-15:sampleend(k),:) = 500;
    else
        final_image_compression(samplestart(k)-2:samplestart(k)+4,:) = 500;
        final_image_compression(sampleend(k)-2:sampleend(k)+4,:) = 500;
    end
end

%% Display
figure;
colormap(gray)
imagesc([-20*.177,20*.177],[0,depth_mm],final_image_compression,[-40,0]) %display image with 40 dB dynamic range
axis image
title(strcat('c=',num2str(c/1000)));%this line will display the speed of sound for you as the figure title

