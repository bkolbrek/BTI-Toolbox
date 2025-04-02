%
% BASS TRANSMISSION INDEX DEMO
%
% This file demonstrates the BTI calculation using the default settings
% that were used in Lara Harris' work:
% 	100 noise iterations
% 	Simple FFT band filtering of the noise
% 	No downsampling
% 	Fixed length noise signals
%
%  --------------------------------------------------- -------------------------------------------
%    This file is part of the Bass Transmission Index (BTI) Toolbox by
%	   Lara Harris and Bjørn Kolbrek
%
%    Copyright (C) 2015-2025 by Lara Harris and Bjørn Kolbrek
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
irfile = 'spen.wav';

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


%% Calculate BTI
bticalc.calculateMTF();

%% Plot the results
Bti = bticalc.getOutputStruct();
plotbti(Bti, 2);
