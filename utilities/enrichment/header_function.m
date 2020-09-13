function [RawResultsPath, DataPath, FigureRoot] =   header_function(DropboxFolder, project) 
    RawResultsPath = [DropboxFolder 'LocalEnrichmentResults\'];
    DataPath = [DropboxFolder 'ProcessedEnrichmentData\' project '\'];
    FigureRoot = [DropboxFolder 'LocalEnrichmentFigures\PipelineOutput\'];
    mkdir(FigureRoot);