function [ScrambledArray1,ScrambledArray2] = scrambleRows(array1,array2)

ScrambledSpotPairs = nan(2,length(array1));
ToFlip = rand(1,length(array2))>0.5;
for i = 1:length(array1)
    DoFlip = ToFlip(i);
    if DoFlip
        ScrambledSpotPairs(1,i) = array2(i);
        ScrambledSpotPairs(2,i) = array1(i);
    else
        ScrambledSpotPairs(1,i) = array1(i);
        ScrambledSpotPairs(2,i) = array2(i);
    end
end
ScrambledArray1 = ScrambledSpotPairs(1,:);
ScrambledArray2 = ScrambledSpotPairs(2,:);

end