function [sig,p,normBool] = sigDiffCalc(a,varargin)

p = inputParser;
addRequired(p,'a',@isvector)
addOptional(p,'b',0,@isvector)
addOptional(p,'paired',false,@logical);

parse(p,a,varargin{:});
a = p.Results.a;
b = p.Results.b;
paired = p.Results.paired;

warning('off','stats:adtest:OutOfRangePLow')
warning('off','stats:adtest:OutOfRangePHigh')

if numel(b)==1
    if adtest(a(:))
        [p,sig] = signrank(a,b);
        normBool = false;
    else
        [sig,p] = ttest(a,b);
        normBool = true;
    end
else
    
    try
        if adtest(cat(1,a(:),b(:)))
            %     disp('non-normal')
            if paired==1 && length(a)==length(b)
                [p,sig] = signrank(a,b);
            else
                [p,sig] = ranksum(a,b);
            end
            normBool = false;
        else
            %     disp('normal')
            if paired==1 && length(a)==length(b)
                [sig,p] = ttest(a,b);
            else
                [sig,p] = ttest2(a,b);
            end
            normBool = true;
        end
    catch
        sig = 0;
        p = NaN;
        normBool = NaN;
    end
end
warning('on','stats:adtest:OutOfRangePLow')
warning('on','stats:adtest:OutOfRangePHigh')
end