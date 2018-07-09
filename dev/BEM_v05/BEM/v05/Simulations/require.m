function []=require(libname,varargin)
global PATH
if(nargin>1)
    Version=varargin{1};
    if(nargin==3)
        bSilent=varargin{2}; 
    else
        bSilent=0;
    end
else
    Version='';
end

if(~isfield(PATH,'STRING'))
    PATH.STRING='';
end
if(~isfield(PATH,libname))
    error(sprintf('%s is not in variable PATH',libname));
else
    folder=getfield(PATH,libname);
    [startIndex, ~, ~, matchStr] = regexp(PATH.STRING, [folder '[a-zA-Z0-9_-]*']);
    if(length(startIndex)>1)
        error('this should not happen');
    else
        if(length(startIndex)==1)
            % Just printing warning of reload
            [~, ~, ~, ~, ~, ~, g]=regexp(matchStr{1},'/');
            oldVersion=g{end};
            if(~bSilent)
                if(~isequal(oldVersion,Version))
                    fprintf('require:\t %s was already loaded with version %s and is now replaced by version %s\n',libname,oldVersion,Version);
                else
                    fprintf('require:\t %s was already loaded with version %s\n',libname,oldVersion);
                end
            end
            oldfolder=[folder oldVersion];
            % removing the path string
            [~, ~, ~, ~, ~, ~, g] = regexp(PATH.STRING, [oldfolder ':']);
            PATH.STRING=strcat(g{:});
            % and the path...
            rmpath(oldfolder);

        end
        folder=[folder Version];
        % now let's do the real loading...
        if(~isdir(folder))
            error(sprintf('Unable to load libray %s version %s',libname,Version));
        else
            addpath(folder);
            if(~bSilent)
%                 fprintf('require:\t addding path %s\n',folder); %%%
%                 UNCOMMENT THIS 
            end
            setfield(PATH,libname,'LOADED');
            PATH.STRING=[PATH.STRING folder ':'];          
        end
    end
end
end

