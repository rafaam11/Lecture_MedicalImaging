%
%	CT project
%   This template is provided to guide you through the CT project
%
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
disp(sprintf('number of rays = %g', nr))
disp(sprintf('number of views = %g', na))

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
figure; imagesc(x, y, phantom');               % NOTE the transpose (') here and the x and y values
colormap('gray'); axis('square'); title('Disk Phantom');
xlabel('Position'); ylabel('Position');

figure; imagesc(r, ang/pi*180, sg1');     % NOTE the transpose (') here and the fact that angle is displayed in degrees
colormap('gray'); title('Sinogram: Disk Phantom');
xlabel('Position (i.e., Rays)'); ylabel('Angle (i.e., Views)');


%- START FROM HERE


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QUESTION F-H - MYSTERY OBJECT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ('Data3_Sinogram2'); %  MISTERY SINOGRAM (sg2 has same format with sg1)
%- START FROM HERE
