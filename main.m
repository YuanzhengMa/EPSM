%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%README
%Adaptive Stolt migration for thermoacoustic imaging
%Here input needs to contain thick and velocity of each layers
%Resolution is improved by interpolation among x-axis and z-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre work
clc
clear all
load bb2
% Define hyper parameters
data = bb; % utrasound data 
fs = 18e6; % sample frequency
tDelay_range = [0e-6,20e-6]; % delay time between pulse transmission and measurement
cc = [2300, 1550]; % vecity vector of n-layer medium [m/s]
thick = [3.7e-2, 5e-2]; % vector of thickness for each layers [m]
sensor_f_band = [0, 10e6]; % [fLow,fHigh] of transducer frequency band
xStep = 300e-6; % distance between two adjoining elements
interpol_method = 'linear'; % Interpolation method for stolt resampling of spectrum, optition 'linear' and 'chirpz'
p = inputParser;% Input parsing object
addParameter(p,'xFftMult',1);% Multiplier for FFT size in x axis.
p.addParameter('yFftMult',1);% Multiplier for FFT size in y axis.
p.addParameter('tFftMult',1);% Multiplier for FFT size in t axis.
p.addParameter('zFftMult',1);% Multiplier for FFT size in z axis.
p.addParameter('upSamp',4);% Multiplier for interpolation FFT size
p.addParameter('hh',1);% Impulse response
p.addParameter('xStart',0);% First x value of scan
p.addParameter('yStart',0);% First x value of scan
validator = @(str) any(strcmp(str,{'linear','chirpz'}));
p.addParameter('interpol','linear',validator);     % Interpolation method
p.addParameter('fc',mean([sensor_f_band(2),sensor_f_band(1)]));           % Center frequency
p.parse;
param = p.Results;                                  % Store results in "param"
% clear p
%% Iteration for imaging with multi-
% allocate Ram
best_bw = imbinarize(255*ones(3232,250));% bw(best image) with lowest entropy
best_im = ones(3232, 250);%best image
% begin iteration , stolt algo needs kx and kz, 
% and omega is represeted by the above two
for delta_t = tDelay_range(1):1e-6:tDelay_range(2)
    [nT, nX] = size(data);% bscan empty
    tPlot = (0:(nT-1))/fs + delta_t;
    xPlot = (0:(nX-1))*xStep;
    % calculate
    %% Calculate dependent variables
    nL = length(thick);                               % Number of layers
    zIF = cumsum([0; thick(:)]);                      % z position of each interface
    dzl = (cc/2)./(sensor_f_band(2)-sensor_f_band(1));% Z resolution for each layer
    zOffset = delta_t*(cc(1)/2);                      % Z offset of meas. window
    %% Set FFT sizes
    nFFTx = param.xFftMult * 2^nextpow2(nX);          % FFT size in x dimension
    nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.hh)-1);    % FFT size time

    %% Create omega and k vectors
    dOmega = ((2*pi*fs)/nFFTt);                             % Omega step size
    omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*dOmega;       % Omega vector, shift to mid point
    omega = ifftshift(omega);                               % Shift as fft output

    kxs = (2*pi)/xStep;                                % Sampling wavenum., x dir.
    kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
    kx = ifftshift(kx); % inver FFT to calculate kx

    %% Calculate bandpass mask for omega
    omegaBandIndex = (omega >= -(2*pi*sensor_f_band(2))) & (omega <= -(2*pi*sensor_f_band(1)));
    nOmegaBand = nnz(omegaBandIndex);% Number of nonzero matrix elements.
    fBand = fs*(nOmegaBand/nFFTt);

    %% Fourier transform along time dimension
    P = fft(data,nFFTt,1);                               % P(omega,x,y)

    %% Matched filtering
    if param.hh ~= 1
        HH = fft(param.hh,nFFTt,1);                         % Impulse resp. spectrum
        P = P.*repmat(conj(HH),[1 nX nY]);            % Matched filt.
    end

    %% Cut out frequency band
    P = P(omegaBandIndex,:,:);
    omegaBand = omega(omegaBandIndex);

    %% Interpolate band to higher resolution (improves accuracy in Stolt interpol.)
    if strcmp(param.interpol,'linear') && param.upSamp ~= 1
        P = fft(ifft(P),param.upSamp*nOmegaBand);
        tmp = (0:(nOmegaBand*param.upSamp-1))*(dOmega/param.upSamp);
        omegaBand = -tmp(end:-1:1) + omegaBand(end);
    end

    %% tDelay compensation
    P = P.*repmat(exp(-1i*omegaBand(:)*delta_t),[1 nX]);

    %% Fourier transform in x and y direction
    Pokxky = fft(P,nFFTx,2);
    clear P

    %% Create "PSM grids" used for extrapolation between interfaces
    [OMEGA_psm,KX_psm] = ndgrid(omegaBand,kx);

    %% Extraplolate wavefield to each interface, and image using Stolt transf.
    im = cell(nL,1);        % Preallocate cell structure for each layer image
    zIm = cell(nL,1);       % Preallocate cell structure for corresponding z axis

    for ii = 1:nL
        disp(['Processing layer ' num2str(ii) ' of ' num2str(nL)]);

        if strcmp(param.interpol,'linear')
            %%  stolt resample with linear interpolation
            nFFTz = param.zFftMult * 2^nextpow2(nT*(sensor_f_band(2)-sensor_f_band(1))/(fs/2));
            kzs = (2*pi)/dzl(ii);                   % Sampling kz
            kz = (2*pi*sensor_f_band(1))/(cc(ii)/2) + (0:(nFFTz-1))*(kzs/nFFTz);    % kz vector
            kz = -kz(end:-1:1);                     % Correspond to neg. omega

            % "Stolt grids" (correspond to low resolution P(kz,kx,ky) wavefield)

            [KZ_st,KX_st] = ndgrid(kz,kx);
            KK_st = (KZ_st.^2 + KX_st.^2);                      % KK^2
            Akzkxky = 1./(1 + (KX_st.^2)./KZ_st.^2);            % Scale factor
            clear KX_st

            KK_st = -sqrt(KK_st .* (KK_st > 0));        % Calc. KK with square root
            Akzkxky(isnan(Akzkxky)) = 0;                % Remove single NaN point.
                                                        % TODO: Check NaN - why?

            % Calc. omega values to be interpolated for
            OMEGA_st = (cc(ii)/2)*KK_st;

            % Interpolate for each (kx,ky)
            Pkzkxky = complex(zeros(nFFTz,nFFTx));
            for jj = 1:nFFTx
                 Pkzkxky(:,jj) = interp1(omegaBand(:),Pokxky(:,jj),...
                    OMEGA_st(:,jj));
            end

            % Values out af range in interpolation are set to NaN. Change to zero.
            nonValid = (isnan(Pkzkxky)) | (OMEGA_st < omegaBand(1)) | ...
                (OMEGA_st > omegaBand(end));
            Pkzkxky(nonValid) = 0;

            % Amplitude scaling (because of variable change from omega to kz)
            Pkzkxky = Pkzkxky.*Akzkxky;

            % If first layer, shift according to measurement z offset
            if ii == 1
                Pkzkxky = Pkzkxky.*exp(1i*KZ_st*zOffset);
            end

            % Inverse transform
            pzxy = ifftn(Pkzkxky);

        else
            %%  Chirp-z appoximation to stolt interpolation
            Ptkxky = ifft(Pokxky,[],1);
            kzc = (2*pi*param.fc)/(cc(ii)/2);

            % Interpolate for each kx
            Pkzkxky = complex(zeros(nOmegaBand,nFFTx));         % Preallocate
            for jj = 1:nFFTx
                % Calc K (scale factor), A (spec. offset) and W (spec. step)
                K = sqrt(1 + (kx(jj)/kzc)^2);
                A = -2*pi*((fLow/fBand)*(1/K-1) + (param.fc/fBand)*(K-1/K));
                W = 2*pi*(1/K)*(fs/fBand)/nFFTt;
                Pkzkxky(:,jj) = qczt(Ptkxky(:,jj),nOmegaBand,W,A);
            end

            % If first layer, shift according to measurement z offset
            if ii == 1
                kz = omegaBand/(cc(ii)/2);
                Pkzkxky = Pkzkxky.*exp(-1i*repmat(kz,1,nFFTx)*zOffset);
            end

            % Inverse transform
            pzxy = ifftn(flip(Pkzkxky,1));

        end

        % Calculate number of image planes in each layer, calc. z axis coordinates
        if ii == 1
            nPlanesZ = ceil((thick(ii)-zOffset)/dzl(ii)) + 1;
            zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zOffset;
        else
            nPlanesZ = ceil(thick(ii)/dzl(ii)) + 1;
            zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zIF(ii);
        end

        % Cut out part of pzxy corresponding to current layer
        im{ii} = pzxy(1:nPlanesZ,1:nX);

        % Migrate to next layer (if not last layer)
        if ii < nL
            KK_psm = OMEGA_psm/(cc(ii)/2);
            
            KZ2_psm = (KK_psm.^2 - KX_psm.^2);
            KZ_psm = sqrt(KZ2_psm .* (KZ2_psm > 0));
            Pokxky = Pokxky .* exp(-1i*KZ_psm*thick(ii));
        end
    end

    %% Assign output
    % im is Focused image
    % xIm is X-axis pixel positions
    xIm = param.xStart + (0:(nX-1))*xStep;
    % zIm is Z-axis pixel positions

    
    %% Plot
    colormap(hot)
    im1_real = abs(im{1});
    im2_real = abs(im{2});
    % normlization 
    im1_real = im1_real./max(max(im1_real));
    im2_real = im2_real./max(max(im2_real));
    im_temp = [im1_real;im2_real];
    im_temp = imadjust(im_temp,[0.1,1],[0,1]);
    
    
    %imagesc([im1_real;im2_real]);
    figure
    imagesc(xIm*1e3,[zIm{1},zIm{2}]*1e3,best_im)
    %set(gca,'CLim',[-10 0])
    axis equal
    ylabel('Z [mm]')
    title('Focused image')
    xlabel('X [mm]')
    %% Find the best img which has the lowest entropy and  show it
       
    im_temp_bw = imbinarize(im_temp);
    if sum(sum(im_temp_bw))<sum(sum(best_bw))
        best_im = im_temp;
        best_bw = im_temp_bw;
    end
    
end
%% Denoise and Plot, here we use deeplearning methods with pretrained network
net = denoisingNetwork('DnCNN');
best_im = denoiseImage(best_im,net);

% Hilbert firlting
im_hilbert = zeros(size(best_im));
Nfir=30;
b=firpm(Nfir,[0.09 0.9],[1 1],'Hilbert');
for j=1:size(best_im,2)      
    y_hilbert=best_im(:,j);
    HOData=conv(y_hilbert,b);
    yh = HOData(Nfir/2+1:size(best_im,1) +Nfir/2);
    if size(yh,1)>size(im_hilbert,1)
        im_hilbert(:,j) = yh(1:size(im_hilbert,1));
    end
    %im_hilbert2(:,j) = yh.^2+y_hilbert.^2;
end
figure

%imagesc([im1_real;im2_real]);
imagesc(xIm*1e3,[zIm{1},zIm{2}]*1e3,best_im)
title('Best Image');
%set(gca,'CLim',[-10 0])
axis equal
ylabel('Z [mm]')
xlabel('X [mm]')