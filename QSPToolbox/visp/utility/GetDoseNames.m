function doseNames = GetDoseNames(doseArray)
 doseNames = cell(1,length(doseArray));
 for theIndex = 1: length(doseArray)
       doseNames{counter}  = get(doseArray{theIndex},'Name');
 end
 for theIndex = 1: length(theDoseIndices)
            counter               = counter + 1;
            theDoseArray          = [theDoseArray, myModelDoseObjects(theDoseIndices(theIndex))];
            intervention(counter) = interventionCounter; 
            names{counter}  = get(myModelDoseObjects(theDoseIndices(theIndex)),'Name');
        end       
end