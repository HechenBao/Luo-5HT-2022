function [output] = getround (inputdata,scale)
    new_length = (round(length(inputdata)./scale)-1).*scale; % avoid exceed index
    output = inputdata(:,1:new_length);
end
