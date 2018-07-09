function [ I ] = whichfile( Files, pattern )
%WHICHFILE returns the index of the files matching the given pattern
    A=regexpi(Files, pattern);
    I=[];
    for i=1:length(A)
        if(~isempty(A{i}))
            I(end+1)=i;
        end
    end
end

