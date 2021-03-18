%% Assignement 5

% must clear all to chage between videos
clear all;
close all;
clc;

%% audio Indicators
WarnWave = [sin(1:0.05:150), sin(1:0.1:150)];
Audio1 = audioplayer(WarnWave, 22050); % done section
WarnWave = [sin(1:0.1:200), sin(1:0.15:200), sin(1:0.2:200)]; % done code 
Audio2 = audioplayer(WarnWave, 22050);
% play(Audio1);
% play(Audio2);

%%  Load Car Vid                                - in use
v = VideoReader('monte_carlo_low.mp4');
video = read(v);
[width height rgbval frames] = size(video);
done = 1
play(Audio1);


%% Load Ski Vid                                 - not in use
% v = VideoReader('ski_drop_low.mp4');
% video = read(v);
% [width height rgbval frames] = size(video);
% done = 1
% play(Audio1);


%% Reshape picture to matrix 
for i = 1:frames
    x = rgb2gray(video(:,:,:,i));
    vidframes(:,i) = double(reshape(x(:,:),width*height,1));
end
done = 2
play(Audio1);


%% Compute SVD
X1 = vidframes(:,1:end-1);
X2 = vidframes(:,2:end);
[U,S,V] = svd(X1,'econ');
s = S;
done = 3
play(Audio1);


%% Extract Necessary Modes
diagS = diag(s)/sum(diag(s))*100;
percent = 0;
i = 1;
while percent < 80
    percent = percent + diagS(i);
    i = i+1;
end
feature = i-1;


U = U(:,1:feature);
S = S(1:feature,1:feature);
V = V(:,1:feature);
done = 4
play(Audio1);


%%  Stild's eig vals are equal to the eig vals of arbitraty matrix A
%   A is the (the koopman opperator)
Stild = U'*X2*V*diag(1./diag(S));
[eigVec, eigVal] = eig(Stild);

Mu = diag(eigVal); % contains eig val of dmd
dt = 1;
omega = log(Mu)/dt; %continous eig val

Phi = U*eigVec; % eig vec of matrix A (the koopman opperator)
% Now we have found eig vec and eig val of A

y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions
t = 1:frames - 1;

omegalater = omega;
omega = omegalater (abs(omega) <= 0.1);
y0 = y0(abs(omegalater) <= 0.1);
Phi = Phi(:,abs(omegalater) <= 0.1);
% Now we have found eig vec and eig val of A

        % in lecture version
u_modes = zeros(length(y0),length(t));
for iter = 1:length(t)
    u_modes(:,iter) = y0.* exp(omega * t(iter));
end
u_dmd = Phi * u_modes;
%   end of DMD alg
done = 5
play(Audio1);


%%
LR = abs(u_dmd);
LRski = LR;
Xsp = X1 - LR;
err = Xsp .* (Xsp < 0);
done = 6
play(Audio1);


%%

    % LR_err = background
LR_err = uint8(LR + err);
    % Xsp_err = foreground
Xsp_err = uint8(Xsp - err);
done = 7
play(Audio1);


%%

    % LR_err = backgrounf 
LR = uint8(LR);
    % Xsp = foreground
Xsp = Xsp_err;
done = 8
play(Audio1);


%% Play Background
close all;

for (i = 1: frames-1)
    image(:,:,i) = reshape(LR(:,i),width,height);
end 

% for (i = 1: frames-1)
%     imshow(image(:,:,i));
%     i
% end 
done = 9
play(Audio1);


%% Before skier forground to make white
% Xsp = X1 - LRski;
% Xsp = (Xsp - min(Xsp))./(max(Xsp) - min(Xsp));
% done = 10
% play(Audio1);


%% Play Foreground
close all;

for (i = 1: frames-1)
    image(:,:,i) = reshape(Xsp(:,i),width,height);
end 
% for (i = 1: frames-1)
%     imshow(image(:,:,i));
%     i
% end
done = 11
play(Audio1);


%% Comparing Snapshots
close all;

snap = 200;
image_og = rgb2gray(video(:,:,:,snap));

% original snap shot
figure(1);
imshow(image_og);

figure(2);
% edited snap shot back 
image = reshape(LR(:,snap),width,height);
imshow(image);

figure(3);
% edited snap shot back 
image = reshape(Xsp(:,snap),width,height);
imshow(image);
done = 12
play(Audio1);


%% Plotting Singular Values
figure(4)
plot((diag(s)/sum(diag(s))*100), 'ob')
set(gca,'Fontsize',18)
title('Single Values')
xlabel('Index of Single Value')
ylabel('Energy of Single Value (%)')
done = 13
play(Audio1);


%% Plotting Omega
figure(5)
plot(real(omegalater), imag(omegalater), 'ob')
set(gca,'Fontsize',18)
title('omega')
xlabel('Real Components of omega')
ylabel('Imaginary Components of omega')
done = 14
play(Audio2);







