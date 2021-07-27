function varargout = opt_args(def_args,varargin)

 

% verify optional arguments

max_opt_args = length(def_args);

num_opt_args = length(varargin);

if num_opt_args > max_opt_args

    ST = dbstack(1);

    error(strcat(ST.name,':TooManyArguments'),'Received %d optional arguments, but expected no more than %d.',num_opt_args,max_opt_args);

end

 

% return argument values

varargout = def_args;

varargout(1:num_opt_args) = varargin;
