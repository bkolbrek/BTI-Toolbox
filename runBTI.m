%
% BASS TRANSMISSION INDEX ANALYSIS
%
% This script has a file open dialog that makes it more convenient for
% normal use when analysing the IR of a loudspeaker. Adjust the parameters
% in the "Parameters" section at the start of the file before running it.
% For rapid testing, using downsampling will speed up the process.
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

% ---------   Parameters   --------- %
numNoiseIterations = 100; % default = 100
minNumModulationCycles = 2; % default = 2
useDownsampling = false; % default = false
apprNewFs = 3000; % approximate new sampling frequency
useConstantNumModulationCycles = false; % default = false;
filterType = ModulationTransferFunctionCalculator.FILTER_BRICKWALL; % default: FILTER_BRICKWALL
% --------- End Parameters --------- %

% remember path
if ~exist('path', 'var')
	path = '';
else
	% safeguard against when the dialog is canceled
	oldpath = path;
end

[irfile, path] = uigetfile('*.wav','Open impulse response', path);

if irfile == 0
	% dialog is canceled, use old path if it exists
	if exist('oldpath', 'var')
		path = oldpath;
	else
		path = '';
	end
	return
end

% Load impulse response
[IR, fs] = audioread([path, irfile]);
IR = IR';
[~, modelName, ~] = fileparts(irfile);
disp('IR file loaded');

% Initialize calculator and load IR
bticalc = ModulationTransferFunctionCalculator();
bticalc.setSystemImpulseResponse(IR, fs, modelName);

% Set calculation parameters
bticalc.setNumNoiseIterations(numNoiseIterations);
bticalc.setMinNumModCycles(minNumModulationCycles);
bticalc.setProcessingSamplingFreq(apprNewFs);
bticalc.setUseDownsampling(useDownsampling);
bticalc.setConstNumModCycles(useConstantNumModulationCycles);
bticalc.setNoiseFilterMethod(filterType);

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

% Calculate BTI
bticalc.calculateMTF();

% Plot the results
Bti = bticalc.getOutputStruct();
plotbti(Bti, 2);
