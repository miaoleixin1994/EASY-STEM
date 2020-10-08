function [ Image ] = ImNorm( Image )
% Just a Image function that brings min to 0 and max to 1;
%   Detailed explanation goes here

Image = (Image - min(Image(:)));

Image = Image/max(Image(:));

end

