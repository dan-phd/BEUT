function avError = relError( varargin )
% Calculates relative error between multiple n-dimensional arrays
% or cells of arrays.
%
% relError( A,B,C,... ) where size(A)=size(B)=size(C) calculates
% the error between all inputs and displays the error for each
% comparison.
%
% avError = relError( A,B,C,... ) outputs a 2D array of size
% nArgs x nArgs comparing each input to each input. Use 'display'
% parameter to suppress the output text.
%
% relError( ..., 'parameter',parameter_value ) where 'parameter'
% can be:
%   'display' to allow printing the result to screen [true]
%   'inf' to set whether or not to allow anomolous infs [false]
%   'NaN' to set whether or not to allow anomolous NaNs [false]
%
% Example
%    A=[1 2 NaN]; B=[1 2 NaN]; C=[2 4 NaN];
%    err = BEUT.relError(A,B,C,'display',false)
%
%    err =
%            A    B    C
%    A -->   0    0   100
%    B -->   0    0   100
%    C -->  50   50    0
%
% So for the error between the first argument and the third
% argument, use err(1,3).
%

%   Author: Daniel Simmons - DansPhD.com
%   Edited: 23/10/2015

% Deal with not enough arguments
nArgs = length(varargin);
if nArgs<2, error('More arguments are needed'); end;

% Parse optional arguments
show=true;
allowNaN = false; allowInf = false;
for a = 1:nArgs
    if ischar(varargin{a})
        
        % Allow user to specify whether an output is displayed
        if strcmpi(varargin{a},'display')
            show=varargin{a+1};
            nArgs=nArgs-2;
            
        elseif strcmpi(varargin{a},'NaN')
            allowNaN=varargin{a+1};
            nArgs=nArgs-2;
            
        elseif strcmpi(varargin{a},'inf')
            allowInf=varargin{a+1};
            nArgs=nArgs-2;
            
        else
            error('Check help for optional parameters supported, not "%s"',varargin{a})
        end
    end
end

% Go though each combination of pairs of input arguments
avError = zeros(nArgs,nArgs);
for k=1:nArgs
    for j=k+1:nArgs
        
        % Deal with cells
        if iscell(varargin{j})
            assert(iscell(varargin{k}),'Each input must be a cell if one is used')
            nk = numel(varargin{k}); nj = numel(varargin{j});
            assert(nj==nk,'Number of elements in all cells must be equal')
            J=varargin{j}; K=varargin{k};
            for cell_num=1:nj
                avError(k,j) = avError(k,j) + compareArray(K{cell_num},J{cell_num});
                avError(j,k) = avError(j,k) + compareArray(J{cell_num},K{cell_num});
            end
        else
            
            avError(k,j) = compareArray(varargin{k}, varargin{j});
            avError(j,k) = compareArray(varargin{j}, varargin{k});
            
        end % if cells
        
        if show
            fprintf('\nRelative error between arrays %i and %i = %.3g %%',...
                k,j,max(avError(j,k),avError(k,j)))
        end
        
    end % j
end % k

% Buffer the display output
if show
    fprintf('\n\n')
end

    % Private function comparing 1 pair of arguments
    function [ av_error ] = compareArray( A, B )
        
        % Deal with matrices with non-matching dimensions
        nA = ndims(A); nB = ndims(B);
        assert(nA==nB,'Matrices must have the same dimensions')
        for i=1:nA
            N = size(A,i); M = size(B,i);
            assert(N==M,'First array has length %i, second has length %i',N,M)
        end
        
        % Translate arguments to 1 dimension
        A = A(:); B = B(:);
        
        % Deal with 0 matrices
        assert(~all(A==0),'Matrices must have values other than 0')
        
        deal_with('inf');
        deal_with('NaN');
        
        
        function deal_with(type)
            
            if strcmp(type,'inf')
                allow_element=allowInf;
                find_fcn = @isinf;
            elseif strcmp(type,'NaN')
                allow_element=allowNaN;
                find_fcn = @isnan;
            end
            
            idx_A = find_fcn(A); idx_B = find_fcn(B);
            
            if any(idx_A) || any(idx_B)
                if all(idx_A == idx_B)
                    A(idx_A) = 0;      % Set elements in both arrays to 0
                    B(idx_B) = 0;      % if both are at the same position
                elseif allow_element
                    A(idx_A) = 0;      % Set inf in both arrays to 0
                    B(idx_A) = 0;      % if anomolous elements are allowed
                    A(idx_B) = 0;
                    B(idx_B) = 0;
                else
                    error(['Array contains anomalous ' type '.'])
                end
            end
        end
        
        
        %% Compute error
        if (all(isreal(A)) && all(isreal(B)))      % If both arrays are real...
%             av_error = 100*mean(abs(A-B))/mean(abs(A));  % Old method
            av_error = 100*( sqrt(sum((A-B).^2)) )/( sqrt(sum(A.^2)) );
        else                                       % If any array has imaginary values...
            av_error = 100*(norm(A-B)/norm(A));
        end
        
    end

end
