function options = paramSet(varargin)
%optionset Create/alter parameters structure.
%   OPTIONS = paramSet('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names. 
%   
%   OPTIONS = paramSet(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   OPTIONS = paramSet(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties. 
%   
%   paramSet with no input arguments displays all property names and their
%   possible values.
%   

default = struct('Isp', [280 300 465 0],...
                 'T', [2000 700 80 0]*1e3,...
				 'Mp', [435 180 32 0]*1e3,...
				 'M', [483 200 35 7]*1e3,...
				 'tb',[129.39 172.1 729.12 inf],...
				 'Sref',pi*(2)^2/4,...
				 'beta0',0,...
				 'g0',9.8);
				 
% Print out possible values of properties.
	if (nargin == 0) && (nargout == 0)
		fprintf('          Isp: [ scalar or row vector (sg) ]\n');
		fprintf('          T:   [ scalar or row vector (N) ]\n');
		fprintf('          tb : [ scalar or row vector (sg) ]\n');
		fprintf('          M  : [ scalar or row vector (kg)  \n');
		fprintf('          Mp : [ scalar or row vector (kg)  \n'); 
		fprintf('\n');
	return;
	end

	Names = [
		'Isp             '
		'T               '
		'tb              '
		'M               '
		'Mp              '
		'Sref            '
		'beta0           '
        'g0              '
		];

	m = size(Names,1);
	names = lower(Names);

	% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
	options = [];
	for j = 1:m
		options.(deblank(Names(j,:))) = default.(deblank(Names(j,:)));%[];
	end
	i = 1;
	while i <= nargin
		arg = varargin{i};
		if ischar(arg) || (isstring(arg) && isscalar(arg)) % arg is an option name
			break;
		end
		if ~isempty(arg)                      % [] is a valid options argument
			if ~isa(arg,'struct')
				error(message('MATLAB:odeset:NoPropNameOrStruct', i));
			end
			for j = 1:m
				if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
					val = arg.(deblank(Names(j,:)));
				else
					val =[];
				end
				if ~isempty(val)
					options.(deblank(Names(j,:))) = val;
				end
			end
		end
		i = i + 1;
	end

% Convert string arguments and options.
	for ii = 1:nargin
		if isstring(varargin{ii}) && isscalar(varargin{ii})
			varargin{ii} = char(varargin{ii});
		end
	end

% A finite state machine to parse name-value pairs.
	if rem(nargin-i+1,2) ~= 0
		error(message('MATLAB:odeset:ArgNameValueMismatch'));
	end
		expectval = 0;                          % start expecting a name, not a value
	while i <= nargin
		arg = varargin{i};
		if ~expectval
			if ~ischar(arg)
				error(message('MATLAB:odeset:NoPropName', i));
			end
			lowArg = lower(arg);
			j = strmatch(lowArg,names);
			if isempty(j)                       % if no matches
				error(message('MATLAB:odeset:InvalidPropName', arg));
			elseif length(j) > 1                % if more than one match
				% Check for any exact matches (in case any names are subsets of others)
				k = strmatch(lowArg,names,'exact');
				if length(k) == 1
					j = k;
				else
					matches = deblank(Names(j(1),:));
					for k = j(2:length(j))'
						matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
					end
					error(message('MATLAB:odeset:AmbiguousPropName',arg,matches));
				end
			end
			expectval = 1;                      % we expect a value next
		else
			options.(deblank(Names(j,:))) = arg;
			expectval = 0;
        end
		i = i + 1;
	end

	if expectval
		error(message('MATLAB:odeset:NoValueForProp', arg));
	end

	options.M0 = fliplr(cumsum(fliplr(options.M)));
	options.t0 = [0,cumsum(options.tb(1:end-1))];
	options.S = 1-options.Mp./options.M;
	options.i = 1:length(options.t0);
	options.r = options.M0./(options.M0-options.Mp);
	options.Rt = 6378e3;
end