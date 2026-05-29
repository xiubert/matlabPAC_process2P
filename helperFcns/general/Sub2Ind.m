function out = Sub2Ind(sz,varargin)
% Sub2Ind: NaN-tolerant 2D linear-index conversion.
%
% Drop-in replacement for MATLAB's built-in sub2ind that lets NaN
% subscripts pass through (the built-in errors on NaN as of R2024a).
% Callers use NaN to mark out-of-bounds neighbors that should remain
% NaN in the output so they can be filtered downstream.
%
% Usage:
%   linIdx = Sub2Ind([m n], iSub, jSub)
%
% Inputs:
%   sz       - size vector [m n] of the target 2D array
%   iSub,jSub- row and column subscripts (any matching shape); NaN
%              entries propagate to NaN in linIdx.
%
% Output:
%   linIdx = iSub + (jSub - 1) * m  (matches sub2ind for finite inputs)
%
% Restrictions:
%   2D only. The per-axis stride below uses sz(k-1) instead of the
%   correct prod(sz(1:k-1)), so ND>2 inputs would yield wrong indices.
%   An explicit check below errors out rather than silently mis-indexing.
%
% See also: sub2ind, ind2sub, anmlFRA2dPrime

  if numel(sz) ~= 2 || numel(varargin) ~= 2
      error('Sub2Ind:NotSupported', ...
          ['Sub2Ind only supports 2D inputs (sz must be length 2 and ' ...
           'exactly 2 subscript arrays must be supplied). Got sz of ' ...
           'length %d and %d subscript inputs.'], numel(sz), numel(varargin));
  end

  out = varargin{1};
  for i = 2:numel(varargin)
      out = out + (varargin{i}-1)*(sz(i-1));
  end
end