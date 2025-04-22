%
% BASS TRANSMISSION INDEX DEMO
%
% This file demonstrates BTI calculation using alternative settings
% to those were used in Lara Harris' work:
% 	Different number of noise iterations
% 	Constant number of modulation periods
% 	Downsampling (requires signal processing toolbox)
% 	Butterworth bandpass filtering of noise (requires DSP System toolbox)
%
% There will be some differences in the results, but usually not anything
% major. Please consult Lara Harris' PhD thesis chapter 2.
%
%  --------------------------------------------------- -------------------------------------------
%    This file is part of the Bass Transmission Index (BTI) Toolbox by
%	   Lara Harris and Bjørn Kolbrek
%
%    Copyright (C) 2025 by Lara Harris and Bjørn Kolbrek
%
%    The BTI Toolbox is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as published by 
%    the Free Software Foundation, version 2.1.
%
%    The BTI Toolbox is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%    FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License along with the
%    BTI Toolbox. If not, see <http://www.gnu.org/licenses/>.
%  --------------------------------------------------- -------------------------------------------
%

% Initialization
clear;
close all;
dataFolder = '.\IRs\';
irfile = 'yam.wav';

% Load impulse response
[IR, fs] = audioread([dataFolder, irfile]);
IR = IR';
[~, modelName, ~] = fileparts(irfile);

% Initialize calculator and load IR
bticalc = ModulationTransferFunctionCalculator();
bticalc.setSystemImpulseResponse(IR, fs, modelName);

% Plot frequency response
[fvec1, frf1] = bticalc.getResponseFunction();
figure(1);
semilogx(fvec1, frf1);
xlim([fvec1(1), fvec1(end)]);
xlabel('Frequency [Hz]')
ylabel('Level [dB]');
grid;
title(modelName);
drawnow;

% Default settings
bticalc.setNumNoiseIterations(100);
bticalc.setMinNumModCycles(2);
bticalc.setUseDownsampling(false);
bticalc.setConstNumModCycles(false);
bticalc.setNoiseFilterMethod(ModulationTransferFunctionCalculator.FILTER_BRICKWALL);


%% Calculate BTI with fewer noise iterations
bticalc.setNumNoiseIterations(10);
bticalc.setUseDownsampling(false);
bticalc.calculateMTF();

% Plot the results
Bti = bticalc.getOutputStruct();
plotbti(Bti, 2);
drawnow;

%% Calculate BTI using a constant number of modulation cycles
bticalc.setNumNoiseIterations(10);
bticalc.setMinNumModCycles(10);
bticalc.setConstNumModCycles(true);
bticalc.calculateMTF();

% Plot the results
Bti = bticalc.getOutputStruct();
plotbti(Bti, 3);
drawnow;

%% Calculate BTI using downsampling
if isOctave()
  pkg load signal;
end
bticalc.setNumNoiseIterations(100);
bticalc.setProcessingSamplingFreq(3000);
bticalc.setUseDownsampling(true);
bticalc.setMinNumModCycles(2);
bticalc.setConstNumModCycles(false);
bticalc.calculateMTF();

% Plot the results
Bti = bticalc.getOutputStruct();
plotbti(Bti, 4);
drawnow;

%% Calculate BTI using Butterworth filtering
% Please note that this may not work in Octave.
bticalc.setSystemImpulseResponse(IR, fs, modelName);
bticalc.setNumNoiseIterations(10);
bticalc.setMinNumModCycles(2);
bticalc.setUseDownsampling(false);
bticalc.setConstNumModCycles(false);
bticalc.setNoiseFilterMethod(ModulationTransferFunctionCalculator.FILTER_BUTTERWORTH);
bticalc.calculateMTF();

% Plot the results
Bti = bticalc.getOutputStruct();
plotbti(Bti, 5);
drawnow;
