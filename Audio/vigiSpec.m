function vigiSpec(varargin)
% vigiSpec - Generates a spectrogram with customizable parameters.
%
% Usage:
%   vigiSpec(w, fs, F, thr, win, step, tR, rgb)
%
% Inputs (in order):
%   w      - Signal (1 x time)
%   fs     - Sampling frequency (Hz)
%   F      - Vector of frequencies for the spectrogram
%   thr    - Scaling parameter (0-1) to adjust noise threshold
%   win    - FFT window size (should be a power of 2)
%   step   - Overlap fraction of the window
%   tR     - Time vector shift range (optional)
%   rgb    - Boolean flag for color freeze mode (0 = normal, 1 = frozen)
%
% Notes:
%   - Default values are optimized for bird song analysis.
%   - If an argument is omitted or empty (`[]`), the default is used.
%   - The function supports automatic range adjustments for colors.

% Set default parameters (Vigi standard for bird song analysis)
defaultArgs = {[], [], 500:5:12e3, 0.6, 256, 180, [], 0};

% Replace missing inputs with default values
optargs = defaultArgs;
optargs(1:nargin) = varargin(1:nargin);
optargs(cellfun(@isempty, optargs)) = defaultArgs(cellfun(@isempty, optargs));
[w, fs, F, thr, win, step, tR, rgb] = optargs{:};

% Compute spectrogram
[~, F, T, P] = spectrogram(w, win, step, F, fs);

% Adjust time range if specified
if ~isempty(tR)
    T = tR;
end

% Convert power spectrum to decibels
S = 10 * log10(P);
S(S == -Inf) = min(S(:)); % Replace -Inf values with the lowest real value

% Set colormap
cmap = jet(256);
cmap(1:8, :) = 0; % Set lowest colors to black

colormap(cmap);

% Display spectrogram
if rgb == 0
    imagesc(T, F, S);
    set(gca, 'YDir', 'normal');
    set(gca, 'CLim', [min(S(:)) + thr * range(S(:)), max(S(:))]); % Adjust color range
else
    % Normalize and scale spectrogram for RGB mode
    map = [min(S(:)) + thr * range(S(:)), 1.1 * max(S(:))];
    S = max(min(S, map(2)), map(1)); % Clip values
    S = ceil((S - map(1)) / range(map) * 255) + 1;
    
    % Convert to RGB image
    Srgb = ind2rgb(S, cmap);
    imagesc(T, F, Srgb);
    set(gca, 'YDir', 'normal');
end
end
