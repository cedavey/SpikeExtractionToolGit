%     out = ternaryOp(cond,a,b,doeval,input)
% If cond is true, return a, else return b
% If you can't evaluate the true or false statement to even run the
% function (e.g. in ternaryOp(isempty(x), 0, x(1)), you can't call the
% function because it needs to evaluate x(1) to get the input parameters,
% but for empty x this will return an error), then input a string & set
% doeval to true, then the evaluation of each statement will only happy if
% the condition is satisfied.
%
% If you would like a or b to simply return the value of the condition
% being evaluated, then set to 'input', e.g. 
%     ternaryOp(cond,'input',b)
% returns 'input' if the condition is true, else returns b.
%
% If the condition is scalar then results a & b can be whatever size, and
% the whole object will be returned. However, if the condition is a vector,
% then either a & b must be vectors, or else match the size of cond.
%
% Inputs:
%     cond - logical condition to check for true/false
%        a - output when condition is true
%        b - output when condition is false
%   doeval - set to true to evaluate a/b only if cond is true/false
%    input - a & b can be function handles, in which case this can be used
%            as an input to the function
function out = ternaryOp(cond,a,b,doeval,input)
   returnInput = 'input';
    cond = logical(cond);
    % assume we have a single test to return all of a or all of b
    if isscalar(cond)
        if nargin<4 || isempty(doeval) || ~doeval
            if cond
               if ischar(a) && strcmpi(a,returnInput)
                  out = cond; 
               else 
                  out = a; 
               end
            else
               if ischar(b) && strcmpi(b,returnInput)
                  out = cond; 
               else 
                  out = b; 
               end
            end
        else
            if nargin<5
                input = [];
            end
            if cond
                out = evalStatement(a,input);
            else
                out = evalStatement(b,input);
            end
        end
        return;
    end    
    
    % assume there's a test per element of cond - if value to set to is a
    % vector only use elmts of the vector where the condition is true/false. 
    if isscalar(b)
       out = repmat(b,size(cond));
    else
       if ~all( size(b) == size(cond) )
          str = 'If using vectors, size of condition must match size of false result\n';
          cprintf('Keyword*', str);
          return;
       end
       out = b;
    end

    if ~isscalar(a)
       % if result a is a vector, and the condition is a vector, sizes must match
       if ~isscalar(cond)
          if ~all( size(a) == size(cond) )
             str = 'If using vectors, size of condition must match size of true result\n';
             cprintf('Keyword*', str);
             return;
          end
          out(cond) = a(cond);
          
       % if result a is a vector, but condition is scalar, simply return the whole vector
       % - we should never get here since this condition is covered above
       else
          if cond
             out = a;
          else
             out = b;
          end
       end
    else
       out(cond) = a;
    end
    if ~isscalar(b)
       out(~cond) = b(~cond);
    else
       out(~cond) = b;
    end
end
        
function out = evalStatement(statement,input)
    % If the condition statement is a function handle, evaluate it
    if isa(statement,'function_handle')
        out = statement(input);
    % If the condition statement is a string, evaluate it
    elseif ischar(input)
        out = eval(statement);
    % Else, output the condition 
    else
        out = statement;
    end
end
