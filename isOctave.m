%isOctave : checks if the code is running on Octave.
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

function io = isOctave()
persistent runningOctave;
if isempty(runningOctave)
    runningOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
end
io = runningOctave;