% -----------------------------------------------------------------------------------
%  Class for calculating the Bass Transmission Index.
%  Based on BTI code by Lara Harris.
%  Refactored into a class by Bjørn Kolbrek.
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

classdef ModulationTransferFunctionCalculator < handle
	properties (Constant)
		FILTER_BRICKWALL = 1;
		FILTER_BUTTERWORTH = 2;
	end

	properties (GetAccess = public, SetAccess = protected)
		noiseMaxAmplitude = 1;
		numNoiseIterations = 100;
		useDefaultRgn = true;
		% useDefaultRgn = true will create a matrix of normalised columns, each column individually normalised.
		% This will make sure we get the same random distribution every time.
		modulationFreqs
		bandLimits
		numBands
		numModFreqs
		noiseNormalized
		noiseBandLimited
		systemName
		systemImpulseResponse
		samplingFreq
		calcSamplingFreq
		origSamplingFreq
		inputEnvelopes
		outputEnvelopes
		tVector
		mtfMatrix
		useDownsampling = false;
		centerImpulseResponse = true;
		lowestModFreq
		noiseLengthInSamples
		noiseFilterMethod
		inputFFTsize
		IRstartIndex
		frequencyResolution = 1.2
		minNumModCycles = 2;
		constNumModCycles = false;
	end

	methods (Access = private)

	end

	methods
		function self = ModulationTransferFunctionCalculator()
			self.calcSamplingFreq = 24000;
			self.init();
		end

		% ---------------------------------------------------------------------
		%
		%   Initialisation methods
		%
		% ---------------------------------------------------------------------
		function init(self)
			self.initBTI();
			self.inputFFTsize = 2^15;
			self.noiseMaxAmplitude = 1;
			self.numNoiseIterations = 100;
			self.numBands = length(self.bandLimits.Centre);
			self.numModFreqs = length(self.modulationFreqs);
			self.systemName = 'New system';
			self.noiseFilterMethod = ModulationTransferFunctionCalculator.FILTER_BRICKWALL;
		end

		function setSystemImpulseResponse(self, systemImpulseResponse, samplingFreq, systemName)
			self.systemImpulseResponse = systemImpulseResponse(:)';
			self.samplingFreq = samplingFreq;
			self.systemName = systemName;
			self.IRstartIndex = self.findIRstartIdx();
		end

		% ---------------------------------------------------------------------
		%
		%   Main calculation method
		%
		% ---------------------------------------------------------------------
		function calculateMTF(self)
			t0 = tic();
			% initialise
			self.initCalculation();
			disp('Initialisation complete:')
			toc(t0)
			% calculate
			if self.constNumModCycles
				% new method, constant same number of modulation cycles
				mfPeriods = self.modulateConstCycles();
			else
				% method from thesis, constant noise length
				mfPeriods = self.modulateConstNoiseLength();
			end
			self.calculateModulationIndex(mfPeriods);
		end

		% This is if you want to plot the response function externally
		function [fvec, frf] = getResponseFunction(self)
			nfft = 2^nextpow2(length(self.systemImpulseResponse));
			fvec = (0:nfft-1) * self.samplingFreq/nfft;
			frf = fft(self.systemImpulseResponse, nfft);
			index = find((fvec >= 10) & (fvec < self.samplingFreq/2));
			fvec = fvec(index);
			frf = 20*log10(abs(frf(index)));
		end

		% the basic info about the results
		function mtfData = getOutputStruct(self)
			mtfData = struct;
			mtfData.SystemName = self.systemName;
			mtfData.Matrix = self.mtfMatrix;
			[~, mtfMeanRounded] = getMtfMean (self);
			mtfData.MeanScore = mtfMeanRounded;
			mtfData.ModulationFreqs = self.modulationFreqs;
			mtfData.CentreFreqs = self.bandLimits.Centre;
			mtfData.BandLimits = [self.bandLimits.Lower, self.bandLimits.Upper];
			mtfData.nNoiseIterations = self.numNoiseIterations;
		end

		% get the mean of the matrix
		function [mtfMean, mtfMeanRounded] = getMtfMean (self)
			%--Overall mean--
			mtfMean = mean( self.mtfMatrix(:) );
			% Decided to roudn to 3dp, so use:
			if isOctave()
				mtfMeanRounded = round(mtfMean*1000)/1000;
			else
				mtfMeanRounded = round(mtfMean, 3);
			end
		end


		% ---------------------------------------------------------------------
		%
		%   Setters
		%
		% ---------------------------------------------------------------------

		function setNumNoiseIterations(self, nIterations)
			self.numNoiseIterations = nIterations;
		end

		function setMinNumModCycles(self, nMin)
			self.minNumModCycles = nMin;
		end

		function setConstNumModCycles(self, constMod)
			self.constNumModCycles = constMod;
		end

		function setProcessingSamplingFreq(self, fs)
			self.calcSamplingFreq = fs;
		end

		function setUseDownsampling(self, useDS)
			self.useDownsampling = useDS;
		end

		function setCenterImpulseResponse(self, cntImp)
			self.centerImpulseResponse = cntImp;
		end

		function setNoiseFilterMethod(self, method)
			self.noiseFilterMethod = method;
		end

	end

	methods (Access = protected)

		function initBTI(self)
			self.modulationFreqs = [0.8, 1.1, 2.2, 4.3, 5.8, 8.5, 11.7];  % As in thesis

			% Define the analysis band frequency limits here, based on thesis by Lara Harris
			nbands = 10;
			inc = 14.8041;  % The difference between adjacent frequencies.

			fL0 = 16.1499;
			fC0 = 22.8041;
			fU0 = 30.9540;

			freqSteps = (0 : nbands-1) * inc;

			f = [fL0, fC0, fU0] + freqSteps';

			self.bandLimits = struct;
			self.bandLimits.Lower = f(:,1);
			self.bandLimits.Centre = f(:,2);
			self.bandLimits.Upper = f(:,3);

		end

		function index = findIRstartIdx(self)
			IR2 = self.systemImpulseResponse.^2;
			[mx, ind] = max(IR2);
			ind0 = find(IR2(1:ind) < mx/10);
			index = ind0(end);
		end

		function downsampleIR(self)
			self.origSamplingFreq = self.samplingFreq;
			ratio = floor(self.samplingFreq / self.calcSamplingFreq);
			self.systemImpulseResponse = decimate(self.systemImpulseResponse, ratio);
			self.samplingFreq = self.samplingFreq / ratio;
			bins = self.frequencyResolution * self.samplingFreq;
			self.inputFFTsize = 2^(nextpow2(bins));
			self.IRstartIndex = self.findIRstartIdx();
		end

		function preconditionIR(self)
			frf = fft(self.systemImpulseResponse, self.inputFFTsize);
			% remove DC offset
			frf(1) = 0;
			% normalise to max level within band limits
			fRes = self.samplingFreq/self.inputFFTsize;
			startIndex = round(self.bandLimits.Lower(1) / fRes) + 1;
			endIndex = round(self.bandLimits.Upper(end) / fRes) + 1;
			magnitudeValidRange = abs(frf(startIndex:endIndex));
			gainFactor = 1/max(magnitudeValidRange);
			frf = gainFactor * frf;
			if self.centerImpulseResponse
				% SHIFT TO MIDDLE

				N = length(self.systemImpulseResponse);
				centreSample = N/2; % what if it's N/2+1??
				[~, m] = max(abs(self.systemImpulseResponse));
				sampleShift = centreSample - m;

				fv = (0:self.inputFFTsize-1)/N; %norm. freq vector
				jexp = -1i * 2*pi.*fv; %Most of the exponential factor
				pf = exp(jexp.*sampleShift); %The full exp factor;

				%-Apply phase shift
				frf = frf.*pf;
			end
			if isOctave()
				self.systemImpulseResponse = ifft(frf, self.inputFFTsize);
			else
				self.systemImpulseResponse = ifft(frf, self.inputFFTsize, 'symmetric');
			end
			self.IRstartIndex = self.findIRstartIdx();
		end

		function initCalculation(self)
			if self.useDownsampling
				self.downsampleIR();
			end
			self.preconditionIR();

			self.lowestModFreq = min(self.modulationFreqs);
			samplesInOnePeriodOfMinModFreq = (1/self.lowestModFreq) * self.samplingFreq;

			self.noiseLengthInSamples = round(samplesInOnePeriodOfMinModFreq * self.minNumModCycles);
			self.tVector = (0:self.noiseLengthInSamples-1) / self.samplingFreq;

			% Preallocate.
			% Each row is different mf. Results will be added for each noise iteration
			% Used cells because it is an array of arrays!
			self.inputEnvelopes = cell(self.numBands, self.numModFreqs);
			self.inputEnvelopes(:) = {zeros(self.noiseLengthInSamples,1)};
			self.outputEnvelopes = self.inputEnvelopes;

			self.createNormalisedNoiseMatrix(self.noiseLengthInSamples);
			self.bandlimitNoise();

			self.mtfMatrix = zeros(self.numBands, self.numModFreqs);
		end

	end

	methods (Access = private)
		function spectrum = getComplexSpectrum(self, nfft)
			spectrum = fft(self.systemImpulseResponse(:), nfft);
		end

		function mfPeriods = modulateConstNoiseLength(self)
			t0 = tic();
			lengthIR = length(self.systemImpulseResponse);
			p = nextpow2(lengthIR + size(self.noiseBandLimited,1) - 1);
			convFftSize = 2^p;  % assumes it's fixed always
			% we only need to find the spectrum of the IR once, as it never changes
			systemSpectrum = self.getComplexSpectrum(convFftSize);

			% --- main noise iteration loop ---
			for noiseCounter = 1 : self.numNoiseIterations
				t1 = tic();
				for mfCounter = 1 : self.numModFreqs
					for bandCounter = 1 : self.numBands
						bandlimitedNoise = self.noiseBandLimited(:,bandCounter,noiseCounter);


						testSignal = self.amplitudeModulate(bandlimitedNoise, self.modulationFreqs(mfCounter));

						outputSignal = self.simulatePlayback(testSignal, systemSpectrum, convFftSize);

						%-Assign outputs---------------------------------------------------
						absOutput =  abs(outputSignal); % OUTPUT SIGNAL- abs
						absInput = abs(testSignal); %INPUT signal mod BL (norm) noise

						% Update data---
						self.outputEnvelopes(bandCounter, mfCounter) = {self.outputEnvelopes{bandCounter, mfCounter} + absOutput};
						self.inputEnvelopes(bandCounter, mfCounter) = {self.inputEnvelopes{bandCounter, mfCounter} + absInput};

					end % EO for all mod freqs
				end % EO for all bands
				display([self.systemName, ': iteration #', num2str(noiseCounter)])
				toc(t1)
			end % EO for all noise iterations
			disp('Modulation complete:')
			toc(t0);

			% info for further processing
			mfPeriods = struct;
			mfPeriods.Secs = 1 ./ self.modulationFreqs;
			mfPeriods.Samples = floor( self.samplingFreq * mfPeriods.Secs);
			mfPeriods.NumSegments = floor(self.noiseLengthInSamples ./ mfPeriods.Samples);
		end

		function mfPeriods = modulateConstCycles(self)
			t0 = tic();

			% info for further processing
			mfPeriods = struct;
			mfPeriods.Secs = 1 ./ self.modulationFreqs;
			mfPeriods.Samples = floor( self.samplingFreq * mfPeriods.Secs);
			mfPeriods.NumSegments(1:self.numModFreqs) = self.minNumModCycles;
			mfPeriods.nfft = mfPeriods.Samples*0;

			lengthIR = length(self.systemImpulseResponse);

			% --- main noise iteration loop ---
			for mfCounter = 1 : self.numModFreqs
				t1 = tic();
				% the size of the fft varies with modulation frequency
				p = nextpow2(lengthIR + mfPeriods.Samples(mfCounter) - 1);
				convFftSize = 2^p;  % assumes it's fixed always
				mfPeriods.nfft(mfCounter) = convFftSize;
				systemSpectrum = self.getComplexSpectrum(convFftSize);
				index = 1:mfPeriods.Samples(mfCounter);
				for noiseCounter = 1 : self.numNoiseIterations
					for bandCounter = 1 : self.numBands
						bandlimitedNoise = self.noiseBandLimited(index,bandCounter,noiseCounter);

						testSignal = self.amplitudeModulate(bandlimitedNoise, self.modulationFreqs(mfCounter));

						outputSignal = self.simulatePlayback(testSignal, systemSpectrum, convFftSize);

						%-Assign outputs---------------------------------------------------
						absOutput =  abs(outputSignal); % OUTPUT SIGNAL- abs
						absInput = abs(testSignal); %INPUT signal mod BL (norm) noise

						% Update data---
						self.outputEnvelopes{bandCounter, mfCounter}(index) = self.outputEnvelopes{bandCounter, mfCounter}(index) + absOutput;
						self.inputEnvelopes{bandCounter, mfCounter}(index) = self.inputEnvelopes{bandCounter, mfCounter}(index) + absInput;


					end % EO for all bands
				end % EO for all noise iterations
				display([self.systemName, ': MF #', num2str(mfCounter)])
				toc(t1)
			end % EO for all mod freqs
			disp('Modulation complete:')
			toc(t0);

		end

		function calculateModulationIndex(self, mfPeriods)
			%{
				Individual iterations to find MTF.
    		* each element in the cell array corresponds to a combination of band and mod freq.

    		* each element contains a 1xN array - the envelopes resulting from
     		  cumulative summing over all noise iterations.

    		* The length of the noise depends on lowest mod freq (fixed to some extent)
    		  and sampling rate - variable but likely limited in practice to a few common
    		  values.
			%}
			t1 = tic();
			% --- Loop for calculating the MTF for all combinations ---
			%FOR each band
			for bandCounter = 1 : self.numBands
				%FOR each mf
				for mfCounter = 1 : self.numModFreqs
					% Extract envelope results.
					outputEnvelopesNow = self.outputEnvelopes{bandCounter, mfCounter};
					inputEnvelopesNow = self.inputEnvelopes{bandCounter, mfCounter};

					nSegments = mfPeriods.NumSegments(mfCounter);
					nSamplesPerSegment = mfPeriods.Samples(mfCounter);  % i.e. samples in one period of current mod frequency

					% Average across all cycles of the ensemble-averaged IP & OP envelopes
					% This is chopping into segments and summing.
					% Reshape and sum along columns. Make sure we use an integer number of cycles,
					% otherwise reshape will not correctly separate cycles.

					outputReshaped = reshape(outputEnvelopesNow(1:nSegments*nSamplesPerSegment), [], nSegments);
					inputReshaped = reshape(inputEnvelopesNow(1:nSegments*nSamplesPerSegment), [], nSegments);

					outputAveragedEnvelope = sum(outputReshaped, 2) / nSegments;
					inputAveragedEnvelope = sum(inputReshaped, 2) / nSegments;


					% Remove offset. ----
					% What is left is the remaining modulation amplitude.
					eOPnow = (outputAveragedEnvelope - min(outputAveragedEnvelope)) ;
					eIPnow = (inputAveragedEnvelope - min(inputAveragedEnvelope)) ;

					%Calculate m from 1 cycle, offsets removed
					modDepthInput = sum(eIPnow);
					modDepthOutput = sum(eOPnow);

					modulationIndex = modDepthOutput/modDepthInput;

					self.mtfMatrix(bandCounter, mfCounter) = modulationIndex;

				end % EO for all modulation frequencies

			end % EO for all frequency bands
			disp('MTF calculation complete:');
			toc(t1);
		end

		function createNormalisedNoiseMatrix(self, numSamples)
			% will create a matrix of normalised columns, each column individually normalised.
			% make sure we get the same random distribution every time
			if self.useDefaultRgn
				rng('default');
			end
			noise = randn(numSamples, self.numNoiseIterations, 'single');
			gainFactor = self.noiseMaxAmplitude ./ max(abs(noise));
			self.noiseNormalized = (noise .* repmat(gainFactor, numSamples, 1));
		end

		function bandlimitNoise(self)
			noiseSize = size(self.noiseNormalized);
			self.noiseBandLimited = zeros(noiseSize(1), self.numBands, self.numNoiseIterations);
			% Filters the noise to the current bandlimits. Needs to be done only
			% once.
			if self.noiseFilterMethod == ModulationTransferFunctionCalculator.FILTER_BUTTERWORTH
				% Butterworth bandpass filter
				bandpassSpecs = fdesign.bandpass;
				bandpassSpecs.Astop1 = 60;
				bandpassSpecs.Astop2 = bandpassSpecs.Astop1;
				passStopFactor = 1.2;

				for bandCounter = 1:self.numBands
					% design filter for the different bands
					bandpassSpecs.Fpass1 = self.bandLimits.Lower(bandCounter) * 2.0 / self.samplingFreq;
					bandpassSpecs.Fstop1 = bandpassSpecs.Fpass1 / passStopFactor;
					bandpassSpecs.Fpass2 = self.bandLimits.Upper(bandCounter) * 2.0 / self.samplingFreq;
					bandpassSpecs.Fstop2 = bandpassSpecs.Fpass2 * passStopFactor;
					bpFilter = design(bandpassSpecs,'butter','matchexactly','passband','SystemObject',true);

					% Filter the noise
					bandlimitedNoise = bpFilter(self.noiseNormalized);
					self.noiseBandLimited(:,bandCounter,:) = bandlimitedNoise(1:self.noiseLengthInSamples,:);
				end
			else
				% Default filtering method: FFT brickwall. Very simple but perhaps
				% not the best.
				%-Noise spectrum-----------------------------------------------
				nfft = self.noiseLengthInSamples;
				%nfft = 2^(nextpow2(self.noiseLengthInSamples));
				noiseResolution = self.samplingFreq/nfft; %max Freq res of noise
				noiseSpectrum = fft(self.noiseNormalized, nfft); %normalised noise

				for bandCounter = 1:self.numBands
					% loop over all bands and filter the noise for all iterations in one go
					%Find band freq limit indexes
					freqIndexLow = round( self.bandLimits.Lower(bandCounter)/noiseResolution) +1; %NOTE +1 because fnz(0) is 0
					freqIndexHigh = round( self.bandLimits.Upper(bandCounter)/noiseResolution) +1;

					%-Creating the bandlimited noise
					noiseSpectrumBandlimited = zeros(size(noiseSpectrum));
					noiseSpectrumBandlimited(freqIndexLow:freqIndexHigh,:) = noiseSpectrum(freqIndexLow:freqIndexHigh,:);
					if isOctave()
						bandlimitedNoise = ifft(noiseSpectrumBandlimited, nfft);
					else
						bandlimitedNoise = ifft(noiseSpectrumBandlimited, nfft, 'symmetric');
					end
					self.noiseBandLimited(:,bandCounter,:) = bandlimitedNoise(1:self.noiseLengthInSamples,:);
				end
			end

		end  % EO function bandlimitnoise

		function modulatedSignal = amplitudeModulate(self, signal, modulationFreq)
			% Create mod function-----------------------------------------
			% < NOTE noise freq resolution must be at least twice the lowest mod freq. to discriminate?>
			modulatingFunction = 0.5 .* (1 + cos((2*pi*modulationFreq)*self.tVector(1:length(signal))));

			modulatedSignal = modulatingFunction(:) .* signal;
		end  % EO function modulatenoise

		function [outputTrimmed] = simulatePlayback(self, signal, systemSpectrum, nfft)
			% simulate passing test signal (noise, bandlimited then modulated) through the loudspeaker:
			% convolve in time domain or multiply in frequency domain.

			N_test = length(signal);
			signalSpectrum = fft(signal(:), nfft);  % Modulated BL noise

			outputSpectrum = signalSpectrum .* systemSpectrum; % Modulated noise x BL IR shifted

			% Back to Time domain
			if isOctave()
				outputSignal = ifft(outputSpectrum);
			else
				outputSignal = ifft(outputSpectrum, 'symmetric');
			end

			%idx1 = round((N_IR / 2))+1;  % MUST do this +1 to get them aligned! Due to centering.
			idx1 = self.IRstartIndex;
			idx2 = idx1 + N_test - 1;
			idxRange = idx1:idx2;  % Should = length of test sig!
			outputTrimmed = outputSignal(idxRange);
		end  % EO function simulateplayback


	end



end
