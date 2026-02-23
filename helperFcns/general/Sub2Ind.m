function out=Sub2Ind(sz,varargin)
  out=varargin{1};
  for i=2:numel(varargin)
      out=out+(varargin{i}-1)*(sz(i-1));
  end
end