# The Bass Transmission Index (BTI) Toolbox

A Matlab toolbox to calculate the Bass Transmission Index for loudspeakers.

(c) Lara Harris and Bjørn Kolbrek

## Support

We do not have the capacity to offer much support for the code or its use, as we are both busy with other projects and work. We will however address actual bugs and errors that are reported. Please read the usage guide, code examples and publications to familiarise yourself with the method and the code. The thesis gives the background and reasons for the choices made in developing the model, and should hopefully answer your questions in that regard. The examples included with the toolbox should be enough to allow you to use the method for your own work. 

If you want to do further research on the BTI (see some suggestions under Future Work below), we would be happy if you made contact with us, although it's unlikely we can offer much in terms of practical involvement. 

## What is the Bass Transmission Index?

The BTI is an objective measure of a loudspeaker's ability to accurately reproduce low frequency content. It was developed by dr. Lara Harris in PhD at the ISVR, University of Southampton, UK, based on previous work by prof. Keith Holland and Philip Newell. The BTI aims to capture how well a loudspeaker reproduces the temporal envelope of a dynamic signal, something that cannot easily be evaluated from the usual frequency response. For the mix engineer, perhaps the biggest problem is how these delays alter the perceived balance between kick drum and bass guitar. A system with poor transient response may cause the mix to be biased towards the kick drum, as the drum transient may not build up to its full level in time, while the resonant bass guitar is enhanced. This is very difficult to correct at a later (mastering) stage, since the two instruments occupy the same part of the spectrum. 

Like the STI, as used in speech intelligibility measurements, the BTI is based on a metric called the Modulation Transfer Function (MTF). MTF-based methods such as these pass an amplitude-modulated signal through the system, varying the rate of the temporal fluctuations across a range of values that are likely to feature in the real-world signals that the system is expected to reproduce. The change in modulation depth will give a measure of how well the system can reproduce these temporal variations.

![Modulation transmission](./Doc/MTF.png)

The BTI evaluates the MTF in 10 frequency bands from 20-160Hz, with modulation frequencies from 0.8-11.7Hz. The flatness of the frequency response magnitude and extension is also part of the metric. The metric produces a value between 0 and 1 for each combination of signal and modulation frequency, which is plotted as a grayscale intensity image. Variation in the horizontal direction indicates variation in the frequency response magnitude, while variation in the vertical direction indicates variation in the faithfulness of envelope transmission. The mean score therefore depends on both low-frequency extension and temporal accuracy.

![BTI](./Doc/BTIplot.png)

Please see dr. Harris' publications for a more in-depth explanation of the background of BTI, including how the method was developed and its parameters selected.

## The Toolbox

This toolbox is based on the code written by Harris for her PhD work. It has been somewhat cleaned up and restructured into a Matlab class by Bjørn Kolbrek. It should be viewed as a reference implementation, and may not be optimal in terms of efficiency or signal processing. 

## License

The BTI Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by 
the Free Software Foundation, version 2.1.

The BTI Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

## Usage

The basic usage is very simple: create an instance of the ModulationTransferFunctionCalculator class, give it an impulse response, and call the calculateMTF() method. This is shown in the BTIdemo.m file. More advanced usage (changing parameters or calculation method) is demonstrated in BTIdemo2.m. 

Also included is the file RunBTI.m, which shows a dialog box allowing you to select any .wav file you want to calculate the BTI of. 

Note that the BTI was developed for use with anechoic IRs. We do not yet have experience with how well the calculation methods in this toolbox work in the presence of reflections etc, but it should not be expected that the results will be accurate. The IR should be long enough to resolve the frequency response down to at least 20Hz with 1-2Hz resolution. Some sort of nearfield, ground plane or outdoor measurements should probably be used if anechoic measurements are not available. 

## How it works

The processing of the impulse response (IR) is done in several stages:
 1. If any downsampling is specified, this is done first. 
 2. The IR is centered.
 3. A matrix of band limited noise is generated.
 4. Looping through all frequency bands and modulation frequencies, the modulated and band limited noise is convolved with the IR.
 5. The modulation transfer function is calculated from the input and output signals.
 6. The results are presented in a grayscale matrix.


## Future work

The BTI has been developed and validated using closed and vented studio monitors measured under anechoic conditions. While the metric shows excellent correlation with the subjective impression of the low frequency performance of these speakers, there are still many use cases that have not been investigated. Some of these are:
 - Room acoustics: Can BTI be used to evaluate speakers in room, or even the low frequency accuracy of the room itself? 
 - Other types of speakers: how well does BTI work with 
   - transmission lines?
   - dipoles?
   - bass horns?
 - Can the calculation be improved further by better signal processing, for instance making it faster?
 - We would also like to have the toolbox translated to other languages like Python and Julia, to make it available to a wider audience.

## Acknowledgements

Thanks to Keith R. Holland for contributing the demo IR files.

## Citation

Cite as: 

L. Harris and B. Kolbrek, 2025, "The BTI Toolbox"

