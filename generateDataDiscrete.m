function data = generateDataDiscrete(opts,optsExtra)

nRuns = opts.nRuns;

if (isfield(optsExtra,'doParforSimple'))
    
    parfor iRun = 1:nRuns
%     for iRun = 1:nRuns
        data(iRun) = InteractingParticlesDiscrete(opts);
    end
    
else
    
    parForStruct = parForLogSetup();
    tempFile = parForStruct.tempFile;
    
    parfor iRun = 1:nRuns
        
        data(iRun) = InteractingParticlesDiscrete(opts);
        
        parForLogUpdateSaveFile(tempFile,nRuns);
        parForLogPrintProgress(parForStruct);
        
    end

    parForLogCleanup(parForStruct.tempDir);

end

if(isfield(opts,'saveHistOnly') && opts.saveHistOnly)
    data = analyseDataDiscrete(data,opts.name,opts.lenValue,false);
    data.opts = opts;
end