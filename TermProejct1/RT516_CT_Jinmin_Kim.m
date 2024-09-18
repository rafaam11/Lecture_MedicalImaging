%   2021 Medical Imaging 
%	CT project 
%   202123008 KIM Jinmin (M.S. Candidate)
%   Department of Robotics Engineering

clear; close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION A-E - KNOWN OBJECT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('Data1_DiskPhantom');
load ('Data2_Sinogram1');

%--------------------------------------------------------
%	Geometry parameters
%--------------------------------------------------------
nr = 128;	dr = 2;		            % number of radial samples and ray spacing
na = nr*2;          	            % number of angular samples
r = dr * ([1:nr]'-(nr+1)/2);	    % radial sample positions
ang = [0:(na-1)]'/na * pi;          % angular sample positions
fprintf('number of rays = %g\n', nr);
fprintf('number of views = %g\n', na);

%--------------------------------------------------------
%	Image parameters: number of pixels, size, etc.
%--------------------------------------------------------
nx = 128; ny = 128;
dx = 2;		                        % 2 mm / pixel
x = dx * ([1:nx]'-(nx+1)/2);
y = -dx * ([1:ny]'-(ny+1)/2);
xx = x(:,ones(1,ny));
yy = y(:,ones(1,ny))';

%--------------------------------------------------------
% Image the Phantom & Sinogram
%--------------------------------------------------------
figure; imagesc(x, y, phantom');       % Figure 1. Disk Phantom              
colormap('gray'); axis('square'); title('Disk Phantom');
xlabel('Position'); ylabel('Position');

figure; imagesc(r, ang/pi*180, sg1');       % Figure 2. Sinogram 1
colormap('gray'); title('Sinogram: Disk Phantom');
xlabel('Position (i.e., Rays)'); ylabel('Angle (i.e., Views)');


%- START FROM HERE
%% Question A
% 0th moment of each projection at 45

theta = 45;                                      % theta = 45�� ���� projection�� �� ����
[~, angle_index] = min(abs(ang-deg2rad(theta))); % 0~180�� 256����� ang���� theta = 45�� �ε����� ã��
projection = sg1(:, angle_index-1);              % ������ ã�� �ε����� �̿��� �� �ε��������� intensity���� sg1���� ã��

figure;     % Figure 3. 0th moment of each projection at 45
project_max=max(projection);                    % scale�� ���� ���
plot(r,projection,'r');                         % theta = 45���� position�� ���� intensity �׷���
axis([r(1),r(nr),0,2*project_max]); title('0th moment of each projection at 45'); % sinogram1�� x���� �����ϰ� ����
xlabel('Position'); ylabel('Intensity');        % label ����

%% Question B
% Simple Backprojection Image
sbp = zeros(nx,ny);       % 128x128�� �� �̹��� ���� sbp (simple backprojection) ����

for i = 1:na              
  temp = zeros(nx);         % 128x128�� �� temp ���� ����         
  for j = 1:nx              % ���� i������ intensity������ ������ temp�� ��� �຤�� ���� ä��
      temp(j, :) = sg1(:, i)'; 
  end
  rotate_angle2 = (i-1)*180/na;    % 180���� ȸ���� ���̹Ƿ� i��° ������ 180/128�� ������
  backproject_theta = imrotate(temp, rotate_angle2,'bicubic','crop');    % temp, rotate_angle�� i������ ��� ����
  sbp = sbp + backproject_theta;    % ���� sbp �̹��� ������ backproject����� ����. i=1~128���� �ݺ�
end

figure;     % Figure 4. Simple Backprojection Image
imagesc(x, y, sbp); colormap('gray'); 
axis('image'); title('Simple Backprojection Image');
xlabel('Position'); ylabel('Position');
  

%% Question C
% Filtered Reconstructed Image

ramlak = abs(linspace(-1,1,nr)');       % -1 ~ 1 ���� ������ 128����Ͽ� ������ ���ؼ� 128x1 Ram-lak Filter�� ������
ramlak = repmat(ramlak, [1 na]);        % Ram-lak ���͸� 256���� Ȯ��
  
sinogramfft = fft(sg1, [], 1);                        % sinogram1�� fft
sinogramfft = fftshift(sinogramfft, 1);             % fftshift�� ���� ��ȯ�� ������ ���ļ��� �߾����� �Ű���
sinogramfre = ifftshift(sinogramfft.*ramlak, 1);    % ramlak���͸� �����Ű�� ��ȯ�� ������ �ٽ� �Ű���
sinogramfilter = real(ifft(sinogramfre, nr, 1));      % ifft�� �ϸ� real�� ���� �Ǽ��θ� ǥ���� 

figure;     % Figure 5. �� = 45�� projection without and with filtering
plot(r, sg1(:,angle_index)./max(sg1(:,angle_index)), 'r-',...
      r, sinogramfilter(:,angle_index)./max(sinogramfilter(:,angle_index)),'b-');
title('�� = 45�� projection without and with filtering'); legend('Without filtering', 'With filtering');
axis([r(1),r(nr),-1,1]); 
xlabel('Position'); ylabel('Projection Value');

figure;     % Figure 6. Filtered Sinogram
imagesc(r, ang/pi*180, sinogramfilter'); colormap('gray'); axis('image');  
title('Filtered Sinogram'); xlabel('Position'); ylabel('Angle');

% Backprojection�� ���� Question B�� ������ ������
fbp = zeros(nx,ny);   
for i = 1:na
  temp = zeros(nx);
  for posy = 1:1:nx
      temp(posy, :) = sinogramfilter(:, i)';
  end
  rotate_angle = ((i-1)*180)/na;
  backproject_theta = imrotate(temp, rotate_angle, 'bicubic','crop'); 
  fbp = fbp + backproject_theta;
end

figure;     % Figure 7. Filtered Reconstructed Image
imagesc(x, y, max(fbp,0)); colormap('gray'); axis('image');  
title('Filtered Reconstructed Image'); xlabel('Position'); ylabel('Position');

%% Question D
% Comparison between 3 curves at y = -7mm line

index = find(y==-7);                    % y = -7�� ���� �ε����� ã��
  line_phantom = phantom(:, index);     % 1. original at y = -7
  line_sbp = sbp(index,:);              % 2. without filtering at y = -7
  line_fbp = fbp(index, :);             % 3. with filtering at y = -7

figure;     % Figure 8. Comparison between 3 curves at y = -7mm line 
plot(y, line_phantom/max(line_phantom),'r-', ...    % max������ ���������μ�, ���� ū���� 1�� �ǵ��� ������
    y, line_sbp/max(line_sbp), 'b-', ...
    y, line_fbp/max(line_fbp),'g-');
legend('Original', 'Simple BP', 'Filtered BP'); title('Comparison between 3 curves at y = -7mm line');
axis([r(1),r(nr),0,1]);
xlabel('Position (mm)'); ylabel('Normalized Value (0 ~ 1)');

%% Question E
% (Downsampling) Filtered Reconstructed Image

Downsample = sg1(:,1:8:end);                % sg1�� �ٿ���ø�. 256���� Column�� 8 �������� 32���� ������.
[dim1, dim2] = size(Downsample);            % dim1 = 128, dim2 = 32
ramlak = abs(linspace(-1,1,nr).');          % 128*1�� Ram-lak ���� ����
ramlak_Down = repmat(ramlak, [1, dim2]);    % Ram-lak ���͸� dim2��ŭ Ȯ�� (128*32)
  
sinogramfft2 = fft(Downsample, [], 1);          % Question C�� ������ ���͸� ���� ����. 
sinogramfft2 = fftshift(sinogramfft2, 1);
sinogramfre2 = ifftshift(sinogramfft2.*ramlak_Down, 1);
sinogramfilter2 = real(ifft(sinogramfre2, nr, 1));
fbp_Down = zeros(nx,ny);

for i = 1:na/8
    temp = zeros(nx);
    for posy = 1:1:nx
        temp(posy, :) = sinogramfilter2(:, i)';
    end
    rotate_angle2 = ((i-1)*180)/(na/8);
    backproject_theta = imrotate(temp, rotate_angle2, 'bicubic','crop'); 
    fbp_Down = fbp_Down + backproject_theta;     
end   

figure;     % Figure 9. (Downsampling) Filtered Reconstructed Image
imagesc(x, y, max(fbp_Down,0)); colormap('gray'); 
axis('image'); title('(Downsampling) Filtered Reconstructed Image');
xlabel('Position'); ylabel('Position');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION F-H - MYSTERY OBJECT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ('Data3_Sinogram2'); %  MISTERY SINOGRAM (sg2 has same format with sg1)
%- START FROM HERE

%% Question F
% Sinogram of unknown object sg2

%--------------------------------------------------------
%	Geometry parameters
%--------------------------------------------------------
nr2 = 367;	dr2 = 1; na2 = 540;         % sg1 => 128*256, sg2 => 367*540      	       
r2 = dr2 * ([1:nr2]'-(nr2+1)/2);	    % radial ���ø�
ang2 = [0:(na2-1)]'/na2 * pi;           % angular ���ø�
fprintf('number of rays = %g\n', nr2);
fprintf('number of views = %g\n', na2);

%--------------------------------------------------------
%	Image parameters: number of pixels, size, etc.
%--------------------------------------------------------
nx2 = 367; ny2 = 367;                   % 367*367�� pixel size�� ����
dx2 = 1;		                        % 1 mm / pixel
x2 = dx2 * ([1:nx2]'-(nx2+1)/2);        % x(Column) ���ø�
y2 = -dx2 * ([1:ny2]'-(ny2+1)/2);       % y(Column) ���ø�
xx2 = x2(:,ones(1,ny2));                % x(Matrix) ���ø�
yy2 = y2(:,ones(1,ny2))';               % y(Matrix) ���ø�

%--------------------------------------------------------
% Image the Phantom & Sinogram
%--------------------------------------------------------
figure;     % Figure 10. Sinogram of unknown object sg2
imagesc(r2, ang2/pi*180, sg2');   
colormap('gray'); title('Sinogram of unknown object sg2');
xlabel('Position'); ylabel('Angle');


%% Question G
% 1. Simple Backprojection Image of sg2
% 2. Filtered Reconstructed Image of sg2

% 1. Question B�� ������ �������, Simple Backprojection Image ����
 sbp2 = zeros(nx2,ny2);
   mid_pt = nx2/2;
   for i = 1:na2
     temp2 = zeros(nx2);
     for posy2 = 1:1:nx2
         temp2(posy2, :) = sg2(:, i)';
     end
     rotate_angle2 = ((i-1)*180)/na2;
     backproject_theta2 = imrotate(temp2, rotate_angle2, 'bicubic','crop'); 
     sbp2 = sbp2 + backproject_theta2;
   end
   
  figure;     % Figure 11. Simple Backprojection Image of sg2
  imagesc(x2, y2, sbp2); colormap('gray'); 
  axis('image'); title('Simple Backprojection Image of sg2');
  xlabel('Position'); ylabel('Position');

  
  
% 2. Question C�� ������ �������, Filtered Reconstructed Image ����
filter2 = abs(linspace(-1,1,nr2).');
filter2 = repmat(filter2, [1 na2]);
sinogramtopad3 = sg2;
sinogramfft3 = fft(sinogramtopad3, [], 1);
sinogramfft3 = fftshift(sinogramfft3, 1);
sinogramfre3 = ifftshift(sinogramfft3.*filter2, 1);
sinogramfilter3 = real(ifft(sinogramfre3, nr2, 1));

fbp2 = zeros(nx2,ny2);
for i = 1:na2
    temp2 = zeros(nx2);
    for posy2 = 1:1:nx2
        temp2(posy2, :) = sinogramfilter3(:, i)';
    end
    rotate_angle2 = ((i-1)*180)/na2;
    backproject_theta2 = imrotate(temp2, rotate_angle2, 'bicubic','crop'); 
    fbp2 = fbp2 + backproject_theta2;
end
    
figure;     % Figure 12. Filtered Reconstructed Image of sg2
imagesc(x2, y2, max(fbp2,0)); colormap('gray'); 
axis('image'); title('Filtered Reconstructed Image of sg2');
xlabel('Position'); ylabel('Position');

%% Question H
% My filter

My_filter = abs((linspace(-1,1,nr2)).^0.5).';         % 367x1 my filter�� ������. ���� Ram-lak���� ������, ������ ����
My_filter = repmat(My_filter, [1 na2]);             % ���͸� 540���� Ȯ��
sinogramfft4 = fft(sg2, [], 1);
sinogramfft4 = fftshift(sinogramfft4, 1);
sinogramfre4 = ifftshift(sinogramfft4.*My_filter, 1);
sinogramfilter4 = real(ifft(sinogramfre4, nr2, 1));

% ���� ����ũ �׸�
figure;     % Figure 13. Mask of Ram-lak filter and My filter
subplot(2,1,1);
plot(r2,filter2); 
title('Ram-lak filter'); axis([-nr2/2,nr2/2,0,1]); 
xlabel('Frequency'); ylabel('Amplitude');
subplot(2,1,2);
plot(r2,My_filter); 
title('My filter'); axis([-nr2/2,nr2/2,0,1]); 
xlabel('Frequency'); ylabel('Amplitude');


% Backprojection process
fbp3 = zeros(nx2,ny2);
for i = 1:na2
    temp2 = zeros(nx2);
    for posy2 = 1:1:nx2
        temp2(posy2, :) = sinogramfilter4(:, i)';
    end
    rotate_angle2 = ((i-1)*180)/na2;
    backproject_theta2 = imrotate(temp2, rotate_angle2, 'bicubic','crop'); 
    fbp3 = fbp3 + backproject_theta2;
end
figure;     % Figure 14. my filtered backprojection
imagesc(x2, y2, max(fbp3,0)); colormap('gray'); 
axis('image'); title('my filtered backprojection');
xlabel('Position'); ylabel('Position');

angle_index = 183;

% Figure 15. Quantitatively evaluate results (at Time domain)
figure;    
plot(r2, sbp2(:,angle_index)./max(sbp2(:,angle_index)),'r-',...
      r2, fbp2(:,angle_index)./max(fbp2(:,angle_index)),'g-',...
      r2, fbp3(:,angle_index)./max(fbp3(:,angle_index)),'b-');
title('Quantitatively evaluate results (at Time domain)'); legend('Simple sg2', 'Filtered sg2', 'My filtered');
axis([r2(1),r2(nr2),-1,1]); 
xlabel('Position'); ylabel('Intensity');

% Figure 16. Quantitatively evaluate results (at Frequency domain)
comp1 = fft(sbp2, [], 1); comp1 = fftshift(comp1, 1);
comp2 = fft(fbp2, [], 1); comp2 = fftshift(comp2, 1);
comp3 = fft(fbp3, [], 1); comp3 = fftshift(comp3, 1);
figure;    
plot(r2, comp1(:,angle_index)./max(comp1(:,angle_index)),'r-',...
      r2, comp2(:,angle_index)./max(comp2(:,angle_index)),'g-',...
      r2, comp3(:,angle_index)./max(comp3(:,angle_index)),'b-');
title('Quantitatively evaluate results (at Frequency domain)'); legend('Simple sg2', 'Filtered sg2', 'My filtered');
axis([-30,30,-1,1]); 
xlabel('Frequency'); ylabel('(Real) Projection Value');


