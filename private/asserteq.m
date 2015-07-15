function asserteq(varargin)
% asserteq(in1, in2, ...)
% Assert all arguments are equal (via "isequal")
% If nargin is 1 then do nothing

% Sam Hallman, 2013

if length(varargin) >= 2
  for i = 2:length(varargin)
    if ~isequal( varargin{1}, varargin{i} )
      error('inputs %d,%d not equal', 1, i);
    end
  end
end
