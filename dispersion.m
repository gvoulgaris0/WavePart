function kh=dispersion(const)
%%
% Solution of the linear dispersion relationship for surface gravity waves.
% The solution satisfies const=kh*tanh(kh) within 1e-06.
% It is obtained using the Newton-Raphson iteration method.
% The input values of const should be zero or positive;
% the routine returns kh = NaN for negative values of const.
%
%% Usage:
%  kh = dispersion(const)
%
%% Input:
% const = vector of omega^2*h/g
% where omega = radian frequency,
%           h = water depth,
%           g = acceleration due to gravity, and
%
%% Output:
%  kh = vector of k*h
%  where k = wavenumber.
%
%% Authors
%  George Voulgaris
%  School of the Earth, Ocean and Environment
%  University of South Carolina, Columbia, SC, USA
%
%% Copyright 2019 George Voulgaris
%
% This file is part of WavePart.
%
% WavePart is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
echo off
%  initialize kh
kh = NaN*ones(size(const));
% solution for zero values of const
m     = find(const==0);
kh(m) = zeros(size(m));
% initial value for positive const
m     = find(const>0);
kh(m) = sqrt(const(m));
m     = find(const>1);
kh(m) = const(m);
% iterative solution for positive const
m     = find(const>0);
while max(abs(const(m) - kh(m).*tanh(kh(m))))>1e-06
    f      = kh(m).*tanh(kh(m)) - const(m);
    fprime = kh(m)./cosh(kh(m)).^2 + tanh(kh(m));
    kh(m)  = kh(m) - f./fprime;
end
end
