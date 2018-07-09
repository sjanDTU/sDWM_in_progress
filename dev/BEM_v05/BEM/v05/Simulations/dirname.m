function [ dirname ] = dirname( file )
%DIRNAME return directory containing file
    A=regexp(file,'/');
    dirname='./';
    if(~isempty(A))
        dirname=file(1:A(end));
    end
end

