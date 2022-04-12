
#configuration
WorkflowPath = os.path.dirname(os.path.realpath(workflow.snakefile))
include: WorkflowPath +'/snakefiles/config.snakefile'
include : WorkflowPath + '/snakefiles/searchgui+pepshaker.snakefile'
#include: WorkflowPath + '/snakefiles/pepGM.snakefile'


#TODO should be possible to simplify the naming scheme a little as name appear in folder & resultsfiles

rule all:
     input:
          #expand(ResultsDir + SampleName + '/Prior{Prior}/'+ HostName + '_' + ReferenceDBName +'_PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',alpha = AlphaRange, beta = BetaRange, Prior = prior),
          ResultsDir +SampleName +'/paramcheck.png'
          #DataDirectory + SampleName + '/' + firstTarget +'_proteinCount.png'


rule CreateFactorGraph:
     input: ResultsDir +'{samplename}/{hostname}_{DBname}_Default_PSM_Report.txt'
     output: ResultsDir +'{samplename}/{hostname}_{DBname}_PepGM_graph.graphml'
     conda: 'envs/graphenv.yml'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName,
          targetTaxa = TargetTaxa,
          firstTarget = firstTarget,
          sourceDB = sourceDB
     log: ResultsDir +'{samplename}/{hostname}_{DBname}.log'
     shell: 'python3 workflow/scripts/CreatePepGMGraph_customJson.py --sourceDB {params.sourceDB} --targetTaxa {params.targetTaxa} --PSM_Report {input} --PeptideMapPath '+ResultsDir+'{params.samplename}/{params.firstTarget}.json --out {output} &>> {log}'

rule RunPepGM:
     input: ResultsDir +'{samplename}/{hostname}_{DBname}_PepGM_graph.graphml'
     output: ResultsDir +'{samplename}/Prior{Prior}/{hostname}_{DBname}_PepGM_Results_a{alpha}_b{beta}_p{Prior}.csv'
     conda: 'envs/graphenv.yml'
     params:
          samplename = SampleName,
          hostname = HostName,
          DBname = ReferenceDBName,
          targetTaxa = TargetTaxa,
          firstTarget = firstTarget
     log : ResultsDir + '{samplename}/Prior{Prior}/log/{hostname}_{DBname}__a{alpha}_b{beta}_pepGM.log'
     shell: 'python3  workflow/scripts/PepGM.py --GraphMLPath {input} --prior {wildcards.Prior} --alpha {wildcards.alpha} --beta {wildcards.beta} --out {output} &>> {log}'

rule BarPlotResults:
     input: ResultsDir +'{samplename}/Prior{Prior}/'+'{hostname}_{DBname}_PepGM_Results_a{alpha}_b{beta}_p{Prior}.csv'
     output:  ResultsDir +'{samplename}/Prior{Prior}/'+'{hostname}_{DBname}_PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',
     conda: 'envs/graphenv.yml'
     params: NumberofResults = TaxaInPlot
     log : ResultsDir + '{samplename}/Prior{Prior}/log/{hostname}_{DBname}__a{alpha}_b{beta}.log'
     shell: 'python3 workflow/scripts/BarPlotResults_customTaxa.py --ResultsFile {input} --NumberofResults {params.NumberofResults} --out {output} &>> {log}'

rule FindBestParameters:
     input: 
          expand(ResultsDir + SampleName + '/Prior{Prior}/'+ HostName + '_' + ReferenceDBName +'_PepGM_Results_a{alpha}_b{beta}_p{Prior}.png',alpha = AlphaRange, beta = BetaRange, Prior = prior)
     output: ResultsDir +'{samplename}/paramcheck.png'
     conda: 'envs/graphenv.yml'
     log : ResultsDir + '{samplename}/paramcheck.log'
     params: Results = ResultsDir+SampleName,
             host = ScientificHostName
     shell: 'python3 workflow/scripts/GridSearchAnalysis.py --host {params.host} --resultsfolder {params.Results} --out {output} &>> {log}'

#rules only needed for debugging/checking results
rule CountProteinsInJson:
     input: ResultsDir +'{samplename}/{taxonname}.json'
     output : ResultsDir + '{samplename}/{taxonname}_proteinCount.png'
     params : NumberofResults = TaxaInProteinCount
     conda: 'envs/graphenv.yml'
     shell: 'python3 workflow/scripts/CountProteins.py --ResultsFile {input} --NumberofResults {params.NumberofResults} --out {output}'