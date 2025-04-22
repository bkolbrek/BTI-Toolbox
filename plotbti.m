function [hFigBti, roundedMatrix] = plotbti(BTI, figno)
% PLOTBTI Plot BTI intensity image from existing results.
% This is based on plotMTF_external. It plots pre-calculated BTI intensity
% images using the struct output from calculatebti or estimatebti.
% Outputs are the figure handle, hFigBti, and numeric results matrix,
% roundedMatrix.
%
% NOTES:
% This is the method used for creating the BTI plots used in publications.
% The function is originally written in Matlab, and the exact looks have been
% a bit difficult to replicate in Octave, but feel free play around with it.
%
% --------------------------------------------------------------------------
%
% Created 24-12-17; Updated 11-07-2019
% Modified to handle Octave by BK 04-03-2025
%
% --------------------------------------------------- ---------------------
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
%  --------------------------------------------------- --------------------
%


%% FONT SETTINGS
fontName = 'Georgia';
% font size is different for Matlab and Octave

%% PLOT BTI- intensity image
hFigBti = figure(figno);

roundingValue = 0.01;
roundingRatio = 1 / roundingValue;
roundedMatrix = round(BTI.Matrix * roundingRatio) / roundingRatio;

nColors = 20;
btiColormap = colormap( gray(nColors) ); % uses n colors from the built-in colormap 'gray'

colormapScalingLimits = [0, 1];
imagesc(roundedMatrix', colormapScalingLimits); %Note for imagesc: cdatamapping is already set to scaled!
axis xy;
axis square
hImgAxes = gca;

colormap(btiColormap)
colorbar('Ylim', colormapScalingLimits);


%% Labels for cf
centreFreqs = BTI.CentreFreqs;

roundedCentreFreqs = num2cell(centreFreqs); % Sure this can be improved. Why do I need rcf??
RCF = cell(1, length(centreFreqs));

for nf = 1 : length(centreFreqs)
	RCF{nf} = num2str(roundedCentreFreqs{nf}, '%-.1f');
end


%% Labels for mf
modulationFreqs = BTI.ModulationFreqs;

roundedModFreqs = num2cell(modulationFreqs);
RMF = cell(1, length(modulationFreqs));

for n = 1 : length(modulationFreqs)
	RMF{n} = num2str(roundedModFreqs{n}, '%-.1f');
end


%% Update plot labels

labelIndexes = 1 : 1 : length(centreFreqs);

set(hFigBti, 'color', [1 1 1])

angleDegrees = 45;
if isOctave()
	fontSize = 18;
	xtickangle(angleDegrees); % Octave-compatible alternative
	set(hImgAxes, 'fontsize', fontSize, 'fontname', fontName, ...
		'tickdir', 'out', ...
		'xTick', labelIndexes,  ...
		'yTick', (1:length(modulationFreqs)), 'yTickLabel', RMF);

	% The code below is required to get the 45 deg rotated x-tick labels in
	% Octave.
	% Erase current tick labels and create new ones
	set(hImgAxes,'XTickLabel',[]);
	xticks = get(hImgAxes,'xtick');
	ylims = get(hImgAxes,'ylim');

	th = text(xticks, repmat( ylims(1)-.05*(ylims(2)-ylims(1)), length(xticks), 1),...
		RCF(labelIndexes),...
		'fontname',fontName,'fontsize',fontSize,...
		'HorizontalAlignment','right','rotation',angleDegrees);

	hXLabel = xlabel('X-Axis Label'); % Set x-axis label
	pos = get(hXLabel, 'position'); % Get current position
	pos(2) = pos(2) - 0.5; % Adjust vertical position (negative moves it downward)
	set(hXLabel, 'position', pos); % Apply new position

	nameString = upper(BTI.SystemName);
	scoreString = num2str(BTI.MeanScore, '%-.3f');

	% Simplified title formatting for Octave, using unicode overbar M
	title([nameString, ': M̅ = ', scoreString]); 

	ylabel('{\it f_m} (Hz)');
	xlabel('{\it f_{cf}} (Hz)');

else
	fontSize = 14;
	hImgAxes.XTickLabelRotation = angleDegrees;
	% don't like yax settings here. should define outside this line?
	set(hImgAxes, 'fontsize',fontSize, 'fontname',fontName,...
		'tickdir','out',...
		'xTick',labelIndexes, 'xTickLabel',RCF(labelIndexes),...
		'yTick',(1 : length(modulationFreqs)), 'yTickLabel', flip(RMF(1:end), 1))

	nameString = upper(BTI.SystemName);
	overbarM = '$$\overline{M}$$';
	scoreString = num2str(BTI.MeanScore, '%-.3f');

	title([nameString, ': ', overbarM, ' = ', scoreString],...
		'Interpreter','latex')
	ylabel('{\it f_m} (Hz)');
	xlabel('{\it f_{cf}} (Hz)') % italics causes the edge of the f to be cut off! hence insert leading space.
end


end %EO functon plotbti.m
